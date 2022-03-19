source("pmwg/sampling.R")

add_info_standard <- function(sampler, prior){
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

fill_samples_standard <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  samples$a_half[, j] <- group_level$a_half
  samples$last_theta_var_inv <- group_level$tvinv
  samples <- fill_samples_base(samples, group_level, proposals, epsilon, j = j, n_pars)
  return(samples)
}

get_startpoints_standard <- function(pmwgs, start_mu, start_var){
  if (is.null(start_mu)) start_mu <- rnorm(pmwgs$n_pars, sd = 1)
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- riwish(pmwgs$n_pars * 3,diag(pmwgs$n_pars))
  start_a_half <- 1 / rgamma(n = pmwgs$n_pars, shape = 2, rate = 1)
  return(list(tmu = start_mu, tvar = start_var, tvinv = MASS::ginv(start_var), a_half = start_a_half))
}

gibbs_step_standard <- function(sampler){
  # Gibbs step for group means, with full covariance matrix estimation
  # tmu = theta_mu, tvar = theta_var
  last <- last_sample_standard(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior
  
  # Here mu is group mean, so we are getting mean and variance
  var_mu <- ginv(sampler$n_subjects * last$tvinv + prior$theta_mu_invar)
  mean_mu <- as.vector(var_mu %*% (last$tvinv %*% apply(last$alpha, 1, sum) +
                                     prior$theta_mu_invar %*% prior$theta_mu_mean))
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  # New sample for mu.
  tmu <- rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
  names(tmu) <- sampler$par_names
  
  # New values for group var
  theta_temp <- last$alpha - tmu
  cov_temp <- (theta_temp) %*% (t(theta_temp))
  if(!is.null(hyper$std_df)){
    B_half <- hyper$std_scale * diag(1, nrow = sampler$n_pars) + cov_temp # nolint
    tvar <- riwish(hyper$std_df + sampler$n_subjects, B_half) # New sample for group variance
    tvinv <- ginv(tvar)
    # Sample new mixing weights.
    a_half <- NULL
  } else{
    B_half <- 2 * hyper$v_half * diag(1 / last$a_half) + cov_temp # nolint
    tvar <- riwish(hyper$v_half + sampler$n_pars - 1 + sampler$n_subjects, B_half) # New sample for group variance
    tvinv <- ginv(tvar)
    
    # Sample new mixing weights.
    a_half <- 1 / rgamma(n = sampler$n_pars,shape = (hyper$v_half + sampler$n_pars) / 2,
                         rate = hyper$v_half * diag(tvinv) + 1/(hyper$A_half^2))
  }
  return(list(tmu = tmu,tvar = tvar,tvinv = tvinv,a_half = a_half,alpha = last$alpha))
}

get_conditionals_standard <- function(s, samples, n_pars, iteration){
  pts2_unwound <- apply(samples$theta_var,3,unwind)
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu,pts2_unwound)
  mu_tilde <- apply(all_samples, 1, mean)
  var_tilde <- var(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[,iteration], unwind(samples$theta_var[,,iteration])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

unwind <- function(var_matrix, ...) {
  y <- t(chol(var_matrix))
  diag(y) <- log(diag(y))
  y[lower.tri(y, diag = TRUE)]
}

last_sample_standard <- function(store) {
  list(
    tmu = store$theta_mu[, store$idx],
    tvar = store$theta_var[, , store$idx],
    alpha = store$alpha[, , store$idx],
    tvinv = store$last_theta_var_inv,
    a_half = store$a_half[, store$idx]
  )
}

sample_store_standard <- function(par_names, subject_ids, iters = 1, stage = "init", ...) {
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  list(
    epsilon = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    origin = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    alpha = array(NA_real_,dim = c(n_pars, n_subjects, iters),dimnames = list(par_names, subject_ids, NULL)),
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    stage = array(stage, iters),
    subj_ll = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    a_half = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL))
  )
}


add_info <- add_info_standard
fill_samples <- fill_samples_standard
get_startpoints <- get_startpoints_standard
gibbs_step <- gibbs_step_standard
get_conditionals <- get_conditionals_standard
sample_store <- sample_store_standard

#Clean up a little bit
rm(add_info_standard)
rm(fill_samples_standard)
rm(get_startpoints_standard)
rm(gibbs_step_standard)
rm(get_conditionals_standard)
rm(sample_store_standard)
