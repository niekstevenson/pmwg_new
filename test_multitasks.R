#Script for testing the recovery of joint models consisting of multiple tasks (can be done in standard or factors or whatever)

rm(list = ls())
source("pmwg/variants/infnt_factor.R")
library(rtdists)

joint_ll <- function(x, data){
  parPreFixs <- gsub("[|].*", "", names(x))
  totalSum <- 0
  for(mod in unique(parPreFixs)){
    currentPars <- x[which(parPreFixs == mod)]
    names(currentPars) <- gsub(".*[|]", "", names(currentPars))
    modelData <- data[[mod]][[1]]
    totalSum <- totalSum + log_likelihood(currentPars, modelData, sample = F)
  }
  return(totalSum)
}

log_likelihood=function(x,data, sample=TRUE) {
  x <- exp(x)
  bPars <- grep("b", names(x))
  bs <- x["A"]+x[bPars][data$condition]
  if (sample) { #for sampling
    out=rLBA(n=nrow(data),A=x["A"],b=bs,t0=x["t0"],mean_v=x[c("v1","v2")],sd_v=c(1,1),distribution="norm",silent=TRUE)
  } else { #for calculating density
    out=dLBA(rt=data$rt,response=data$resp,A=x["A"],b=bs,t0=x["t0"],mean_v=x[c("v1","v2")],sd_v=c(1,1),distribution="norm",silent=TRUE)
    bad=(out<1e-10)|(!is.finite(out))
    out[bad]=1e-10
    out=sum(log(out))
  }
  out
}


n.trials <- 100      #number trials per subject per conditions
n.subj <- 50 #number of subjects
n.cond <- 2
n.exp <- 3

allparameters <- numeric()
alldata <- list()


for(i in 1:n.exp){
  names=c("subject","rt","resp","condition") #names of columns
  data <- data.frame(matrix(NA, ncol = length(names), nrow = (n.trials*n.subj*n.cond))) #empty data frame
  names(data) <- names
  data$condition <- rep(1:n.cond,times = n.trials) #filling in condition
  data$subject <- rep(1:n.subj, each = n.trials*n.cond) #filling in subjects
  alldata[[i]] <- data
  parameter.names <- c(c("A","v1","v2","t0"), c(paste0("b", 1:n.cond)))
  parameters <- c((0.1 + .02*i)*(-1)^i + c(0.4, 1.2, .7, -2.4), seq(0.1, 0.5, length.out = n.cond))
  names(parameters) <- paste0("Mod", i, "|", parameter.names)
  allparameters <- c(allparameters, parameters)
}

fillCorMatrix <- function(matrix, name1, name2, value){
  matrix[grep(name1, rownames(matrix)), grep(name2, colnames(matrix))] <- value
  matrix[grep(name2, rownames(matrix)), grep(name1, colnames(matrix))] <- value
  return(matrix)
}

# adjustCorMatrix <- function(matrix, adj){
#   #Was probably a cleverer way to do this
#   rowNames <- rownames(matrix)
#   for(i in 1:nrow(matrix)){
#     rowPreFix <- substr(rowNames[i], 0, 4)
#     matrix[i,!grepl(rowPreFix, colnames(matrix))] <- .3*matrix[i,!grepl(rowPreFix, colnames(matrix))] 
#   }
#   return(matrix)
# }
# 
# exp(allparameters)
# n.parameters=length(allparameters)
# corMat <- matrix(0, nrow = n.parameters, ncol = n.parameters)
# rownames(corMat) <- colnames(corMat) <- names(allparameters)
# corMat[grep("A", rownames(corMat)),] <- .3 
# corMat[,grep("A", colnames(corMat))] <- .3 
# 
# corMat <- fillCorMatrix(corMat, "v", "b", .25)
# corMat <- fillCorMatrix(corMat, "v", "v", .4)
# corMat <- fillCorMatrix(corMat, "b", "b", .4)
# corMat <- fillCorMatrix(corMat, "t0", "v", -.1)
# corMat <- fillCorMatrix(corMat, "t0", "b", -.25)
# corMat <- adjustCorMatrix(corMat, .33)


# 
# vars = abs(allparameters)/10 #off diagonal correlations are done as absolute/10
# 
# ###std dev correlation on diagonal - you might think this should be corr = 1, but it's actually the standard deviation 
# diag(corMat)=sqrt(vars)
# corMat <- apply(burned$samples$theta_var, 1:2, mean) #sdcor2cov(corMat)
n_factors <- 2
n_pars <- length(allparameters)

num_eff <- round(n_pars/4) + sample(1:(round(n_pars/2)), n_factors)
lambda <- matrix(0, n_pars, n_factors)
for(k in 1:n_factors){
  num_eff_k <- sample(1:n_pars, num_eff[k])
  lambda[num_eff_k, k] = rnorm(num_eff[k], mean = 0, sd = .3)
}


sigma <- 1/rgamma(n_pars, 20, 4) * diag(n_pars)
covMat <- lambda %*% t(lambda) + sigma
subj_random_effects <- t(mvtnorm::rmvnorm(n.subj,mean=allparameters,sigma=covMat))
exp(subj_random_effects)
parPreFixs <- gsub("[|].*", "", rownames(subj_random_effects))
df <- data.frame(subject = 1:n.subj)
for(i in 1:n.exp){
  modIdx <- which(parPreFixs == paste0("Mod", i))
  task_random_effects <- subj_random_effects[modIdx,]
  rownames(task_random_effects) <- gsub(".*[|]", "", rownames(task_random_effects))
  data <- alldata[[i]]
  for (j in 1:n.subj){
    tmp<- log_likelihood(task_random_effects[,j],sample=TRUE,data=data[data$subject==j,])
    data$rt[data$subject==j]=tmp$rt
    data$resp[data$subject==j]=tmp$response
  }
  df[paste0("Mod", i)] <- I(list(split(data, f = data$subject)))
}


pars <- rownames(subj_random_effects)
# Create the Particle Metropolis within Gibbs sampler object ------------------
sampler <- pmwgs(
  data = df,
  pars = pars,
  ll_func = joint_ll,
)
# start the sampler ---------------------------------------------------------
sampler <- init(sampler, n_cores = 12) # i don't use any start points here

# Sample! -------------------------------------------------------------------
burned <- run_stage(sampler, stage = "burn",iter = 500, particles = 100, n_cores =10, pstar = .6, components = c(6, 12, 18))
adapted <- run_stage(burned, stage = "adapt", iter = 5000, particles = 100, n_cores = 10, pstar = .6)
sampled <- run_stage(adapted, stage = "sample", iter = 1000, particles = 100, n_cores = 10, pstar = .6)

corrplot::corrplot(cov2cor(covMat))
cov_recov <- apply(burned$samples$theta_var[,,400:500], 1:2, mean)
corrplot::corrplot(cov2cor(cov_recov))

corrplot::corrplot(cov2cor(covMat))

# N factors heuristic -----------------------------------------------------
library(factor.switching)
lambda_recov <- burned$samples$theta_lambda
n_burn <- 400
criterion <- .05
f_keep <- colMeans(abs(apply(lambda_recov[,,-(1:n_burn)], 1:2, mean))) > criterion
lambda_recov <- lambda_recov[,f_keep, -(1:n_burn)]


# Post process ------------------------------------------------------------

lambda_reordered <- matrix(aperm(lambda_recov,perm = c(2,1,3)), prod(dim(lambda_recov)[1:2]), dim(lambda_recov)[3], byrow = F)
names_all <- expand.grid(1:n_factors, 1:n_pars)
rownames(lambda_reordered) <- paste0("LambdaV", names_all[,2], "_", names_all[,1])
lambda_reordered <- t(lambda_reordered)

lambda_post <- rsp_exact( lambda_mcmc = lambda_reordered,
                          maxIter = 100,
                          threshold = 1e-6,
                          verbose=TRUE )

# Reorders lambda matrix
reshuffle_loadings <- function(lambda, part = .2){
  # Redo the sign
  if(length(dim(lambda)) == 3){
    mean_lambda <- apply(lambda, 1:2, mean)
    sign <- colMeans(mean_lambda) < 0
    lambda[,sign,] <- -1 * lambda[,sign,]
  } else{
    sign <- colMeans(lambda) < 0
    lambda[,sign] <- -1 * lambda[,sign]
    
  }
  
  # Redo the columns order
  n_factors <- ncol(lambda)
  mean_n_loadings <- numeric(n_factors)
  n_pars <- 1
  if(length(dim(lambda)) == 3){    
    mean_lambda <- apply(lambda, 1:2, mean)
    order_idx <- order(colMeans(abs(mean_lambda)), decreasing = T)
    lambda <- lambda[,order_idx,]
  } else{
    order_idx <- order(colMeans((abs(lambda))), decreasing = T)
    lambda <- lambda[,order_idx]
  }
  
  return(lambda)
}

# Making sure they're in the same order and same sign
match_loadings <- function(lambda_true, lambda_recov){
  find_index <- function(values, indexes){
    j <- 0
    notFound <- T
    while(notFound){
      j <- j + 1
      index <- order(values)[j]
      if(!(index %in% indexes)) notFound <- F
    }
    return(list(index = index, value = values[index]))
  }
  # Redo the sign
  dim3 <- F
  if(length(dim(lambda_recov)) == 3){
    dim3 <- T
    lambda_full <- lambda_recov
    lambda_recov <- apply(lambda_full, 1:2, mean)
  }
  
  lambda_true <- reshuffle_loadings(lambda_true)
  k_true <- ncol(lambda_true)
  order <- numeric(k_true)
  sign <- numeric(k_true)
  lambda_copy <- lambda_recov
  for(k in 1:k_true){
    diff_pos <- colMeans(abs(lambda_copy - lambda_true[,k]))
    diff_neg <- colMeans(abs(-1 * lambda_copy - lambda_true[,k]))
    min_pos <- find_index(diff_pos, order)
    min_neg <- find_index(diff_neg, order)
    if(min_pos$value < min_neg$value){
      order[k] <- min_pos$index
      sign[k] <- 1
    } else{
      order[k] <- min_neg$index
      sign[k] <- -1
    }
  }
  if(dim3){
    lambda_recov <- lambda_full[,order,]
    for(k in 1:k){
      lambda_recov[,k,] <- lambda_recov[,k,] * sign[k]
    }
  } else{
    lambda_recov <- lambda_recov[,order]
    for(k in 1:k){
      lambda_recov[,k] <- lambda_recov[,k] * sign[k]
    }
  }
  return(list(lambda_true = lambda_true, lambda_recov = lambda_recov))
}

lambdas <- match_loadings(lambda_true = lambda, lambda_recov = lambda_post$lambda_hat)


# Plotting ----------------------------------------------------------------

corrplot::corrplot(lambdas$lambda_true, is.corr = F, cl.pos = "n")
corrplot::corrplot(lambdas$lambda_recov, is.corr = F, cl.pos = "n")

n_plot <- 18
corrplot::corrplot(cov2cor(covMat)[1:n_plot, 1:n_plot])
corrplot::corrplot(cov2cor(apply(burned$samples$theta_var, 1:2, mean)[1:n_plot, 1:n_plot]))

N_total <- burned$samples$idx - n_burn
lambda_reorder_array <- array(NA_real_, dim = c(n_pars, n_factors, N_total))
for(i in 1:N_total){
  lambda_reorder_array[,,i] <- matrix(lambda_post$lambda_reordered_mcmc[i,], nrow = n_pars, ncol = n_factors, byrow = T)
}

lambdas <- match_loadings(lambda_true = lambda, lambda_recov = lambda_reorder_array)

library(bayesplot)
library(ggplot2)
par_names <- paste0("V", 1:50)
n_pars_used <- 18
for(k in 1:ncol(lambdas$lambda_recov)){
  true <- lambdas$lambda_true[1:n_pars_used,k]
  par_names_k <- paste0(par_names, "_", k)
  mcmc_out <- t(lambdas$lambda_recov[1:n_pars_used,k,])
  colnames(mcmc_out) <- par_names_k[1:n_pars_used]
  print(mcmc_recover_hist(true = true, mcmc_out))
}



