rm(list = ls())
source("pmwg/sampling_multiGroup.R")
library(rtdists)

fillCorMatrix <- function(matrix, name1, name2, value){
  matrix[grep(name1, rownames(matrix)), grep(name2, colnames(matrix))] <- value
  matrix[grep(name2, rownames(matrix)), grep(name1, colnames(matrix))] <- value
  return(matrix)
}

log_likelihood=function(x,data, sample=F) {
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


n.trials <- 125      #number trials per subject per conditions
n.subj <- 25 #number of subjects
n.cond <- 2
n.group <- 2


df <- data.frame()
for(i in 1:n.group){
  names=c("subject","rt","resp","condition") #names of columns
  data <- data.frame(matrix(NA, ncol = length(names), nrow = (n.trials*n.subj*n.cond))) #empty data frame
  names(data) <- names
  data$condition <- rep(1:n.cond,times = n.trials) #filling in condition
  data$subject <- rep(1:n.subj, each = n.trials*n.cond) #filling in subjects
  parameters <- c((0.1 + 0.05*i)*(-1)^i + c(c(0.4, 1.2, .7, -2.4), seq(0.1, 0.5, length.out = n.cond)))

  n.parameters <- length(parameters)
  names(parameters) <- c(c("A","v1","v2","t0"), c(paste0("b", 1:n.cond)))
  print(exp(parameters))
  corMat <- matrix(0, nrow = n.parameters, ncol = n.parameters)
  rownames(corMat) <- colnames(corMat) <- names(parameters)
  corMat <- fillCorMatrix(corMat, "v", "b", .25)
  corMat <- fillCorMatrix(corMat, "v", "v", .4)
  corMat <- fillCorMatrix(corMat, "b", "b", .4)
  corMat <- fillCorMatrix(corMat, "t0", "v", -.1)
  corMat <- fillCorMatrix(corMat, "t0", "b", -.25)
  corMat[grep("A", rownames(corMat)),] <- .3 
  corMat[,grep("A", colnames(corMat))] <- .3 
  n.parameters=length(parameters)
  rownames(corMat) <- colnames(corMat) <- names(parameters)
  vars = abs(parameters)/10
  ###std dev correlation on diagonal - you might think this should be corr = 1, but it's actually the standard deviation 
  diag(corMat)=sqrt(vars)
  corMat <- sdcor2cov(corMat)
  subj_random_effects <- t(mvtnorm::rmvnorm(n.subj,mean=parameters,sigma=corMat))
  for (j in 1:n.subj){
    tmp<- log_likelihood(subj_random_effects[,j],sample=TRUE,data=data[data$subject==j,])
    data$rt[data$subject==j]=tmp$rt
    data$resp[data$subject==j]=tmp$response
  }
  data$group <- i
  data$subject <- data$subject+(i-1)*n.subj
  df <- rbind(df, data)
}

source("pmwg/sampling_multiGroup.R")
pars <- rownames(subj_random_effects)
priors <- list(
  theta_mu_mean = rep(0, length(pars)),
  theta_mu_var = diag(rep(1, length(pars)))
)
sampler <- pmwgs(
  data = df,
  pars = pars,
  prior = priors,
  ll_func = log_likelihood
)
sampler <- init(sampler, n_cores = 15) # i don't use any start points here
# Sample! -------------------------------------------------------------------
burned <- run_stage(burned, stage = "burn",iter = 500, particles = 100, n_cores = 15, pstar = .7)
adapted <- run_stage(burned, stage = "adapt", iter = 100, particles = 100, n_cores = 1, pstar =.4)
sampled <- run_stage(adapted, stage = "sample", iter = 1000, particles = 100, n_cores = 15, pstar = .6)
