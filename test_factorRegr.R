#Script for testing the recovery of factor pmwg s

rm(list = ls())
source("pmwg/sampling_factorRegression.R")

ll <- function(pars, data){
  n_pars <- length(pars)
  return(sum(dmvnorm(as.matrix(data[,2:(1 + n_pars)]), pars, diag(1, nrow = n_pars), log = T)))
}

n_pars <- 15 # j = 1... p; # p with j = 1... p 
n_subjects <-100; # n
n_trials <- 200      #number trials per subject per conditions
n_factors <- 2
n_preds <- 3

parNames <- paste0("V", 1:n_pars)

lambda <- matrix(rep(c(.5,0,.8,.4,-.3), n_factors*3), nrow = n_pars, ncol = n_factors, byrow = TRUE)

sigma <- diag(rep(c(0.15, 0.25, 0.3, 0.2, 0.1), 3))
x <- rmvnorm(n_subjects, mean = c(0, 0, 0), sigma = diag(1, n_preds))
colnames(x) <- paste0("P", 1:n_preds)

beta <- matrix(c(-.4, -.15, .5, .5, .3,
                 .8), nrow = n_preds, byrow = T)

eta <- x %*% beta + rmvnorm(n_subjects, sigma = diag(rep(1, n_factors)))
epsilon <- rmvnorm(n_subjects, sigma = sigma)
y <- sweep(eta %*% t(lambda) + epsilon, 2, -seq(-2,2, length.out = n_pars))

all_data <- data.frame()
for(j in 1:n_subjects){
  sub_data <- data.frame(subject = rep(j, n_trials))
  sub_data[,2:(n_pars + 1)] <- data.frame(rmvnorm(n_trials, y[j,], diag(1, nrow = n_pars)))
  all_data <- rbind(all_data, sub_data)
}

# Create the Particle Metropolis within Gibbs sampler object ------------------
source("pmwg/sampling_factorRegression.R")
sampler <- pmwgs(
  data = all_data,
  pars = parNames,
  ll_func = ll,
  n_factors = n_factors,
  covariates = x
)
# start the sampler ---------------------------------------------------------
sampler <- init(sampler, n_cores = 8) # i don't use any start points here

# Sample! -------------------------------------------------------------------
burned <- run_stage(sampler, stage = "burn",iter = 1000, particles = 150, n_cores =1, pstar = .6)
save(burned, file = paste0("samples/factorRegr_", n_factors, "F_", n_subjects, "S_", n_pars, "P_FactRecovery.RData"))
adapted <- run_stage(burned, stage = "adapt",iter = 1500, particles = 150, n_cores =8, pstar = .6)
sampled <- run_stage(adapted, stage = "sample",iter = 1500, particles = 150, n_cores =8, pstar = .6)
save(sampled, file = paste0("samples/factorRegr_", n_factors, "F_", n_subjects, "S_", n_pars, "P_FactRecovery.RData"))


