rm(list = ls())
source("pmwg/variants/nested.R")

n.trials <- 100
n.subj <- 10
n.groups <- 10
n.pars <- 4


fillCorMatrix <- function(matrix, name1, name2, value){
  matrix[grep(name1, rownames(matrix)), grep(name2, colnames(matrix))] <- value
  matrix[grep(name2, rownames(matrix)), grep(name1, colnames(matrix))] <- value
  return(matrix)
}

ll <- function(pars, data){
  return(sum(dmvnorm(as.matrix(data[,3:(2 + n.pars)]), pars, diag(1, nrow = n.pars), log = T)))
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

pop_args <- list(
  variant = "pmwg/variants/factor.R",
  n_factors = 2
)

group_args <- list(
  variant = "pmwg/variants/diag.R"
)

sampler <- pmwgs(
  data = all_data,
  pars = pars,
  ll_func = ll,
  pop_lvl = pop_args,
  group_lvl = group_args
)

sampler <- init(sampler, n_cores = 10) # i don't use any start points here
# Sample! -------------------------------------------------------------------
source("pmwg/variants/nested.R")
burned <- run_stage(sampler, stage = "burn",iter = 250, particles = 100, n_cores = 10, pstar = .7)
adapted <- run_stage(burned, stage = "adapt", iter = 1000, particles = 100, n_cores = 15, pstar =.7, min_unique = 100)
sampled <- run_stage(adapted, stage = "sample", iter = 500, particles = 100, n_cores = 15, pstar = .7)
