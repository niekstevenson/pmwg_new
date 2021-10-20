rm(list = ls())
source("pmwg/sampling_nested.R")
sourceCpp("pmwg/utilityFunctions.cpp")


n.trials <- 100
n.subj <- 20
n.groups <- 10
n.pars <- 4


fillCorMatrix <- function(matrix, name1, name2, value){
  matrix[grep(name1, rownames(matrix)), grep(name2, colnames(matrix))] <- value
  matrix[grep(name2, rownames(matrix)), grep(name1, colnames(matrix))] <- value
  return(matrix)
}

ll <- function(pars, data){
  return(sum(dmvnrm_arma_fast(as.matrix(data[,3:(2 + n.pars)]), pars, diag(1, nrow = n.pars), log = T)))
}

top_par_mu <- c(1:n.pars)
names(top_par_mu) <- paste0("V", 1:n.pars)

#Fill in some values
top_corMat <- matrix(c(.4, .25, .1, -.3,
                      .25, .6, -.1, .18,
                      .1, -.1, .5, .45,
                      -.3, .18, .45, .3), 
                 nrow = n.pars, ncol = n.pars)
rownames(top_corMat) <- colnames(top_corMat) <- names(top_par_mu)
top_corMat <- sdcor2cov(top_corMat)

all_data <- data.frame(subject = rep(1:(n.subj*n.groups), each = n.trials), group = rep(1:n.groups, each = n.subj*n.trials),
                       v1 = numeric(n.subj*n.groups*n.trials), v2 = numeric(n.subj*n.groups*n.trials), 
                       v3 = numeric(n.subj*n.groups*n.trials), v4 = numeric(n.subj*n.groups*n.trials))
group_data <- data.frame(group = 1:n.groups, v1 = numeric(n.groups), v2 = numeric(n.groups),
                         v3 = numeric(n.groups), v4 = numeric(n.groups))
group_covs <- array(0, dim = c(n.pars, n.pars, n.groups))

seqs <- c(-.5, .5, length.out = n.subj)
for(i in 1:n.groups){
  group_mean <- rmvnorm(1, top_par_mu, top_corMat)
  A <- matrix(runif(n.pars^2)*2-1, ncol=n.pars)
  group_cov <- A %*% t(A)
  for(j in 1:n.subj){
    sub_mean <- rmvnorm(1, group_mean, group_cov)
    sub_data <- rmvnorm(n.trials, sub_mean, diag(1, nrow = n.pars))
    all_data[all_data$subject == (j + n.subj*(i-1)),3:(2+n.pars)] <- sub_data
  }
  group_data[i, 2:(1+n.pars)] <- group_mean
  group_covs[,,i] <- group_cov
}

colMeans(all_data)

pars <- names(top_par_mu)

priors <- list(
  theta_mu_mean = rep(0, length(pars)),
  theta_mu_var = diag(rep(1, length(pars)))
)

sampler <- pmwgs(
  data = all_data,
  pars = pars,
  prior = priors,
  ll_func = ll
)

source("pmwg/sampling_nested.R")
sampler <- init(sampler, n_cores = 15) # i don't use any start points here
# Sample! -------------------------------------------------------------------
burned <- run_stage(sampler, stage = "burn",iter = 500, particles = 100, n_cores = 15, pstar = .7)
adapted <- run_stage(burned, stage = "adapt", iter = 100, particles = 100, n_cores = 1, pstar =.4)
sampled <- run_stage(adapted, stage = "sample", iter = 1000, particles = 100, n_cores = 15, pstar = .6)
