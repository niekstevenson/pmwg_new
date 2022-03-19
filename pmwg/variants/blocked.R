source("pmwg/sampling.R")
source("pmwg/variants/standard.R")
library(magic)

add_info_blocked <- function(sampler, prior, ...){
  # blocked specific attributes
  par_groups <- list(...)$par_groups
  sampler$par_groups <- par_groups
  sampler$n_par_groups <- length(unique(par_groups))
  # Checking and default priors
  if (is.null(prior)) {
    prior <- list(theta_mu_mean = rep(0, sampler$n_pars), theta_mu_var = diag(rep(1, sampler$n_pars)))
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- ginv(prior$theta_mu_var) #Inverse of the matrix
  
  #Hyper parameters
  attr(sampler, "v_half") <- 2
  attr(sampler, "A_half") <- 1
  sampler$prior <- prior
  return(sampler)
}

get_startpoints_blocked <- function(pmwgs, start_mu, start_var){
  if (is.null(start_mu)) start_mu <- rnorm(pmwgs$n_pars, sd = 1)
  # If no starting point for group var just sample some
  if (is.null(start_var)) {
    start_var <- matrix(nrow = 0, ncol = 0)
    for(i in 1:pmwgs$n_par_groups){
      # Check how many parameters are in the current group
      n_pars_group <- sum(pmwgs$par_groups == i) 
      # Make a subblock of start variances for those parameters
      start_var <- adiag(start_var, riwish(n_pars_group * 3,diag(n_pars_group)))  
    }
  }
  start_a_half <- 1 / rgamma(n = pmwgs$n_pars, shape = 2, rate = 1)
  return(list(tmu = start_mu, tvar = start_var, tvinv = MASS::ginv(start_var), a_half = start_a_half))
}

gibbs_step_blocked <- function(sampler){
  # Gibbs step for group means, with full covariance matrix estimation
  # tmu = theta_mu, tvar = theta_var
  tmu_out <- numeric(sampler$n_pars)
  a_half_out <- numeric(sampler$n_pars)
  alpha_out <- sampler$samples$alpha[,,sampler$samples$idx]
  tvar_out <- matrix(0, nrow = sampler$n_pars, ncol = sampler$n_pars)
  tvinv_out <- matrix(0, nrow = sampler$n_pars, ncol = sampler$n_pars)
  idx <- 0
  for(group in 1:sampler$n_par_groups){
    group_idx <- sampler$par_groups == group
    last <- last_sample_blocked(sampler$samples, group_idx)
    hyper <- attributes(sampler)
    prior <- sampler$prior
    
    n_pars <- sum(group_idx)
    
    # Here mu is group mean, so we are getting mean and variance
    var_mu <- ginv(sampler$n_subjects * last$tvinv + prior$theta_mu_invar[group_idx, group_idx])
    mean_mu <- as.vector(var_mu %*% (last$tvinv %*% apply(last$alpha, 1, sum) +
                                       prior$theta_mu_invar[group_idx, group_idx] %*% prior$theta_mu_mean[group_idx]))
    chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
    # New sample for mu.
    tmu <- rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
    names(tmu) <- sampler$par_names[group_idx]
    
    # New values for group var
    theta_temp <- last$alpha - tmu
    cov_temp <- (theta_temp) %*% (t(theta_temp))
    B_half <- 2 * hyper$v_half * diag(1 / last$a_half) + cov_temp # nolint
    tvar <- riwish(hyper$v_half + n_pars - 1 + sampler$n_subjects, B_half) # New sample for group variance
    tvinv <- ginv(tvar)
    
    # Sample new mixing weights.
    a_half <- 1 / rgamma(n = n_pars,shape = (hyper$v_half + n_pars) / 2,
                         rate = hyper$v_half * diag(tvinv) + hyper$A_half)
    tmu_out[(idx+1):(idx+n_pars)] <- tmu
    a_half_out[(idx+1):(idx+n_pars)] <- a_half
    tvar_out[(idx+1):(idx+n_pars), (idx+1):(idx+n_pars)] <- tvar
    tvinv_out[(idx+1):(idx+n_pars), (idx+1):(idx+n_pars)] <- tvinv
    idx <- idx + n_pars
  }
  names(tmu_out) <- sampler$par_names
  return(list(tmu = tmu_out,tvar = tvar_out,tvinv = tvinv_out,a_half = a_half_out,alpha = alpha_out))
}


get_conditionals_blocked <- function(s, samples, n_pars, iteration){
  pts2_unwound <- apply(samples$theta_var,3,unwind)
  index <- rowMeans(pts2_unwound == 0) == 0 #remove all 0 elements
  pts2_unwound <- pts2_unwound[index]
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu,pts2_unwound)
  mu_tilde <- apply(all_samples, 1, mean)
  var_tilde <- var(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[,iteration], unwind(samples$theta_var[,,iteration])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

last_sample_blocked <- function(store, group_idx) {
  list(
    tmu = store$theta_mu[group_idx, store$idx],
    tvar = store$theta_var[group_idx, group_idx,store$idx],
    alpha = store$alpha[group_idx,, store$idx],
    tvinv = store$last_theta_var_inv[group_idx,group_idx],
    a_half = store$a_half[group_idx,store$idx]
  )
}

add_info <- add_info_blocked
get_startpoints <- get_startpoints_blocked
gibbs_step <- gibbs_step_blocked
get_conditionals <- get_conditionals_blocked

#Clean up a little bit
rm(add_info_blocked)
rm(get_startpoints_blocked)
rm(gibbs_step_blocked)
rm(get_conditionals_blocked)

