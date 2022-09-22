#Script for testing the recovery of factor pmwg s

rm(list = ls())
source("pmwg/variants/infnt_factor.R")

ll <- function(pars, data){
  n_pars <- length(pars)
  return(sum(dnorm(t(as.matrix(data[,2:(1 + n_pars)])), pars, 1, log = T)))
}

# Generate data -----------------------------------------------------------
n_pars <- 50 # j = 1... p; # p with j = 1... p 
n_subjects <-50 # n
n_factors <- 3
n_trials <- 150

num_eff <- round(n_pars/4) + sample(1:(round(n_pars/2)), n_factors)
lambda <- matrix(0, n_pars, n_factors)
for(k in 1:n_factors){
  num_eff_k <- sample(1:n_pars, num_eff[k])
  lambda[num_eff_k, k] = rnorm(num_eff[k], mean = 0, sd = .3)
}


sigma <- 1/rgamma(n_pars, 20, 4) * diag(n_pars)
covMat <- lambda %*% t(lambda) + sigma
y <- rmvnorm(n_subjects, mean = seq(-3,3, length.out = n_pars), covMat)


all_data <- data.frame()
for(j in 1:n_subjects){
  sub_data <- data.frame(subject = rep(j, n_trials))
  sub_data[,2:(n_pars + 1)] <- data.frame(rmvnorm(n_trials, y[j,], diag(1, nrow = n_pars)))
  all_data <- rbind(all_data, sub_data)
}

parNames <- paste0("V", 1:n_pars)
##Factor 
# Create the Particle Metropolis within Gibbs sampler object ------------------
source("pmwg/variants/infnt_factor.R")

sampler <- pmwgs(
  data = all_data,
  pars = parNames,
  ll_func = ll
)

n_cores <- 10

# start the sampler ---------------------------------------------------------
sampler <- init(sampler, n_cores = n_cores) # i don't use any start points here

# Sample! -------------------------------------------------------------------
burned <- run_stage(sampler, stage = "burn",iter = 1500, particles = 75, n_cores = n_cores , pstar = .6)
save(burned, file = paste0("samples/factor_", n_factors, "F_", n_subjects, "S_", n_pars, "P_FactRecovery_burned.RData"))
adapted <- run_stage(burned, stage = "adapt",iter = 1500, particles = 100, n_cores = n_cores, pstar = .6, min_unique = 200, n_cores_conditional = n_cores)
sampled <- run_stage(adapted, stage = "sample",iter = 100, particles = 100, n_cores = n_cores, pstar = .6, n_cores_conditional = n_cores)
sampled <- run_stage(sampled, stage = "sample",iter = 700, particles = 100, n_cores = n_cores, pstar = .6, n_cores_conditional = n_cores)

# N factors heuristic -----------------------------------------------------
library(factor.switching)
lambda_recov <- burned$samples$theta_lambda
n_burn <- 500
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

n_plot <- 50
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
n_pars_used <- 20
for(k in 1:ncol(lambdas$lambda_recov)){
  true <- lambdas$lambda_true[1:n_pars_used,k]
  par_names_k <- paste0(par_names, "_", k)
  mcmc_out <- t(lambdas$lambda_recov[1:n_pars_used,k,])
  colnames(mcmc_out) <- par_names_k[1:n_pars_used]
  print(mcmc_recover_hist(true = true, mcmc_out))
}

