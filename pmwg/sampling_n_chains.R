### This is the pmwg sampler for standard puproses. It includes a multivariate normal with full covariance matrix on the group level. 

library(MASS) ## For matrix inverse.
library(MCMCpack) #For inverse wishart
library(lme4)
library(parallel)
library(mvtnorm) ## For the multivariate normal.
library(condMVNorm)

source("pmwg/utils_n_chains.R")
source("pmwg/messaging.R")

pmwgs <- function(data, pars, ll_func, prior = NULL, hyper = NULL, n_chains) {
  ###Gets/sets priors, creates pmwgs object and stores essentials
  # Descriptives
  n_pars <- length(pars)
  subjects <- unique(data$subject)
  n_subjects <- length(subjects)
  
  # Hyperparameters
  v_half <- 2 # hyperparameter on Σ prior (Half-t degrees of freedom)
  A_half <- 1 # hyperparameter on Σ prior (Half-t scale) #nolint
  
  # Storage for the samples.
  samples <- replicate(n_chains, sample_store(pars, subjects), simplify = F)
  # Checking and default priors
  if (is.null(prior)) {
    prior <- list(theta_mu_mean = rep(0, n_pars), theta_mu_var = diag(rep(1, n_pars)))
  }
  # Things I save rather than re-compute inside the loops.
  if(is.matrix(prior$theta_mu_var)){
    prior$theta_mu_invar <- ginv(prior$theta_mu_var) #Inverse of the matrix
  } else{
    prior$theta_mu_invar <- 1/prior$theta_mu_var #Inverse of the matrix
  }
  sampler <- list(
    data = data,
    par_names = pars,
    n_pars = n_pars,
    n_subjects = n_subjects,
    subjects = subjects,
    prior = prior,
    ll_func = ll_func,
    samples = samples,
    n_chains = n_chains,
    init = FALSE
  )
  #Hyper parameters
  attr(sampler, "v_half") <- v_half
  attr(sampler, "A_half") <- A_half
  if(!is.null(hyper$std_shape)){
    attr(sampler, "std_shape") <- hyper$std_shape
    attr(sampler, "std_rate") <- hyper$std_rate
  }
  if(!is.null(hyper$std_df)){
    attr(sampler, "std_df") <- hyper$std_df
    attr(sampler, "std_scale") <- hyper$std_scale
  }
  class(sampler) <- "pmwgs"
  sampler
}

init <- function(pmwgs, start_mu = NULL, start_var = NULL,
                 display_progress = TRUE, particles = 1000, n_cores = 1, epsilon = NULL) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  if (is.null(start_mu)) start_mu <- rnorm(pmwgs$n_pars, sd = 1)
  # If no starting point for group var just sample some
  if (is.null(start_var)) {
    #If prior on covariances is a vector, assume diagonal only, otherwise assume full cvs structure
    if(is.matrix(pmwgs$prior$theta_mu_var)){
      start_var <- riwish(pmwgs$n_pars * 3,diag(pmwgs$n_pars))
    } else{
      start_var <- diag(1/rgamma(pmwgs$n_pars, 10, 5)) #But stupid maybe as startpoint
    }
  }
  
  # Sample the mixture variables' initial values.
  a_half <- 1 / rgamma(n = pmwgs$n_pars, shape = 2, rate = 1)
  # Create and fill initial random effects for each subject
  likelihoods <- array(NA_real_, dim = c(pmwgs$n_subjects))
  if(n_cores > 1){
    proposals <- mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals,start_mu = start_mu, 
                          start_var = start_var, n_particles = particles, pmwgs = pmwgs, 
                          n_chains = pmwgs$n_chains, mc.cores = n_cores)
  } else{
    proposals <- lapply(X=1:pmwgs$n_subjects,FUN=start_proposals,start_mu = start_mu, 
                        start_var = start_var, n_particles = particles, 
                        n_chains = pmwgs$n_chains, pmwgs = pmwgs)
  }
  proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects, pmwgs$n_chains))
  
  tvinv = ginv(start_var)
  group_start <- rep(list(list(tmu = start_mu, tvar = start_var, tvinv = tvinv, a_half = a_half)), pmwgs$n_chains)
  
  pmwgs$samples <- lapply(X = 1:pmwgs$n_chains, fill_samples, 
                          samples = pmwgs$samples,
                          group_level = group_start,
                          proposals = proposals,
                          epsilon = rep(set_epsilon(pmwgs$n_pars, epsilon), pmwgs$n_subjects))
  
  pmwgs$init <- TRUE
  return(pmwgs)
}

fill_samples <- function(chain_idx, samples, group_level, proposals, epsilon, j = 1){

  samples <- samples[[chain_idx]]
  proposals <- proposals[,,chain_idx]
  group_level <- group_level[[chain_idx]]
  samples$theta_mu[, j] <- group_level$tmu
  samples$theta_var[, , j] <- group_level$tvar
  samples$last_theta_var_inv <- group_level$tvinv
  samples$a_half[, j] <- group_level$a_half
  
  samples$alpha[, , j] <- proposals[1:length(group_level$tmu),]
  samples$subj_ll[, j] <- proposals[length(group_level$tmu) + 1,]
  samples$idx <- j
  samples$epsilon[,j] <- epsilon
  samples$origin[,j] <- proposals[length(group_level$tmu) + 2,]
  return(samples)
}


start_proposals <- function(s, start_mu, start_var, n_particles, pmwgs, n_chains){
  #Draw the first start point
  proposals <- particle_draws(n_particles*n_chains, start_mu, start_var)
  colnames(proposals) <- rownames(pmwgs$samples$theta_mu) # preserve par names
  lw <- apply(proposals,1,pmwgs$ll_func,data = pmwgs$data[pmwgs$data$subject == pmwgs$subjects[s], ])
  output <- vector("list", length = n_chains)
  for(i in 1:n_chains){
    lw_tmp <- lw[(n_particles*(i-1)+1):(i*n_particles)]
    weight <- exp(lw_tmp - max(lw_tmp))
    idx <- sample(x = n_particles, size = 1, prob = weight)
    output[[i]] <- list(proposal = proposals[idx,], ll = lw[idx], origin = 2)
  }
  return(output)
}

gibbs_step_diag<- function(samples, sampler, prior, hyper){
  # Gibbs step for diagonal only
  # Get single iter versions, tmu = theta_mu, tvar = theta_var
  last <- last_sample(sampler$samples)
  
  prior$theta_mu_invar <- diag(prior$theta_mu_invar)
  last$tvinv <- diag(last$tvinv)
  
  #Mu
  var_mu = 1.0 / (sampler$n_subjects * last$tvinv + prior$theta_mu_invar)
  mean_mu = var_mu * ((apply(last$alpha, 1, sum) * last$tvinv + prior$theta_mu_mean * prior$theta_mu_invar))
  tmu <- rnorm(sampler$n_pars, mean_mu, sd = sqrt(var_mu))
  names(tmu) <- sampler$par_names
  
  if(!is.null(hyper$std_shape)){
    # InvGamma alternative (probably inferior) prior
    shape = hyper$std_shape + sampler$n_subjects / 2
    rate = hyper$std_rate + rowSums( (last$alpha-tmu)^2 ) / 2
    tvinv = rgamma(n=sampler$n_pars, shape=shape, rate=rate)
    tvar = 1/tvinv
    a_half <- NULL
  } else {
    tvinv = rgamma(n=sampler$n_pars, shape=hyper$v_half/2 + sampler$n_subjects/2, rate=hyper$v_half/last$a_half + 
                     rowSums( (last$alpha-tmu)^2 ) / 2)
    tvar = 1/tvinv
    #Contrary to standard pmwg I use shape, rate for IG()
    a_half <- 1 / rgamma(n = sampler$n_pars, shape = (hyper$v_half + sampler$n_pars) / 2,
                         rate = hyper$v_half * tvinv + 1/hyper$A_half)
  }
  return(list(tmu = tmu, tvar = diag(tvar), tvinv = diag(tvinv), a_half = a_half, alpha = last$alpha))
}

gibbs_step_full <- function(samples, sampler, prior, hyper){
  # Gibbs step for group means, with full covariance matrix estimation
  # tmu = theta_mu, tvar = theta_var
  last <- last_sample(samples)
  
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
                         rate = hyper$v_half * diag(tvinv) + hyper$A_half)
  }
  return(list(tmu = tmu,tvar = tvar,tvinv = tvinv,a_half = a_half,alpha = last$alpha))
}

new_particle <- function (s, data, num_particles, parameters, eff_mu = NULL, 
                          eff_var = NULL, mix_proportion = c(0.5, 0.5, 0), 
                          likelihood_func = NULL, epsilon = NULL, subjects, n_chains = 1) 
{
  particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
  cumuNumbers <- cumsum(particle_numbers)
  eff_mu <- eff_mu[, s]
  eff_var <- eff_var[, , s]
  proposals <- matrix(0, nrow = num_particles * n_chains, ncol = length(parameters[[1]]$tmu))
  for(i in 1:n_chains){
    pars_tmp <- parameters[[i]]
    mu <- pars_tmp$tmu
    var <- pars_tmp$tvar
    subj_mu <- pars_tmp$alpha[, s]
    
    pop_particles <- particle_draws(particle_numbers[1], mu, 
                                    var)
    ind_particles <- particle_draws(particle_numbers[2], subj_mu, 
                                    var * epsilon[s]^2)
    pop_particles[1, ] <- subj_mu
    if(mix_proportion[3] == 0){
      eff_particles <- NULL
    } else{
      eff_particles <- particle_draws(particle_numbers[3], eff_mu, eff_var)
    }
    proposals[(num_particles*(i-1) +1):(num_particles*i),] <- rbind(pop_particles, ind_particles, eff_particles)
    
  }
  #Crucial to do this bit in parallel
  colnames(proposals) <- names(mu)
  lw <- apply(proposals, 1, likelihood_func, data = data[data$subject==subjects[s],])
  
  output <- vector("list", length = n_chains)
  for(i in 1:n_chains){
    pars_tmp <- parameters[[i]]
    mu <- pars_tmp$tmu
    var <- pars_tmp$tvar
    subj_mu <- pars_tmp$alpha[, s]
    
    lw_tmp <- lw[(num_particles*(i-1) +1):(num_particles*i)]
    proposals_tmp <- proposals[(num_particles*(i-1) +1):(num_particles*i),]
    lp <- dmvnorm(x = proposals_tmp, mean = mu, sigma = var, 
                  log = TRUE)
    prop_density <- dmvnorm(x = proposals_tmp, mean = subj_mu, 
                            sigma = var * (epsilon[s]^2))
    if (mix_proportion[3] == 0) {
      eff_density <- 0
    }
    else {
      eff_density <- dmvnorm(x = proposals_tmp, mean = eff_mu, sigma = eff_var)
    }
    lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] * 
                                               prop_density) + (mix_proportion[3] * eff_density))
    l <- lw_tmp + lp - lm
    weights <- exp(l - max(l))
    idx <- sample(x = num_particles, size = 1, prob = weights)
    origin <- min(which(idx <= cumuNumbers))
    output[[i]] <- list(proposal = proposals[idx, ], ll = lw_tmp[idx], origin = origin)
  }
  return(output)
}

get_conditionals_diag <- function(s, samples, n_pars, iteration){
  pts2_unwound <- apply(samples$theta_var,3,diag)
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu,pts2_unwound)
  mu_tilde <- apply(all_samples, 1, mean)
  var_tilde <- var(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[,iteration], diag(samples$theta_var[,,iteration])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

get_conditionals_full <- function(s, samples, n_pars, iteration){
  pts2_unwound <- apply(samples$theta_var,3,unwind)
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu,pts2_unwound)
  mu_tilde <- apply(all_samples, 1, mean)
  var_tilde <- var(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[,iteration], unwind(samples$theta_var[,,iteration])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

run_stage <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 1000,
                      display_progress = TRUE,
                      n_cores = 1,
                      n_unique = ifelse(stage == "adapt", 20, NA),
                      epsilon = NULL,
                      pstar = NULL,
                      mix = NULL,
                      pdist_update_n = ifelse(stage == "sample", 500, NA),
                      epsilon_upper_bound = 2) {
  # Set defaults for NULL values
  mix <- set_mix(stage, mix)
  # Set necessary local variables
  .n_unique <- n_unique
  # Set stable (fixed) new_sample argument for this run
  n_pars <- length(pmwgs$par_names)
  
  if(is.matrix(pmwgs$prior$theta_mu_var)){
    gibbs_step <- gibbs_step_full
    get_conditionals <- get_conditionals_full
  } else{
    gibbs_step <- gibbs_step_diag
    get_conditionals <- get_conditionals_diag
  }
  
  # Display stage to screen
  msgs <- list(
    burn = "Phase 1: Burn in\n",
    adapt = "Phase 2: Adaptation\n",
    sample = "Phase 3: Sampling\n"
  )
  cat(msgs[[stage]])
  
  alphaStar=-qnorm(pstar/2) #Idk about this one
  n0=round(5/(pstar*(1-pstar))) #Also not questioning this math for now
  if(is.null(epsilon)){
    epsilon <- pmwgs$samples$epsilon[,ncol(pmwgs$samples$epsilon)]
  }
  if(length(epsilon) == 1){
    epsilon <- rep(epsilon, pmwgs$n_subjects)
  }
  epsilon <- rep(.3, pmwgs$n_subjects)
  
  n_chains <- pmwgs$n_chains
  # Build new sample storage
  pmwgs$samples <- lapply(pmwgs$samples, extend_sampler, pmwgs, iter, stage)
  
  # create progress bar
  if (display_progress) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples[[1]]$idx
  
  eff_mu <- NULL
  eff_var <- NULL
  data <- pmwgs$data
  subjects <- pmwgs$subjects
  prior <- pmwgs$prior
  hyper <- attributes(pmwgs)
  # Main iteration loop
  for (i in 1:iter) {
    if (display_progress) {
      accRate <- mean(accept_rate(pmwgs))
      update_progress_bar(pb, i, extra = accRate)
    }
    # Create/update efficient proposal distribution if we are in sampling phase.
    if(stage == "sample" & (i %% pdist_update_n == 0 || i == 1)){
      test_samples <- extract_samples(pmwgs, stage = c("adapt", "sample"))
      iteration <- dim(test_samples$theta_mu)[2]
      conditionals=lapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, 
                          n_pars, iteration)
      conditionals <- simplify2array(conditionals)
      eff_mu <- do.call(cbind, conditionals[1,])
      eff_var <- abind(conditionals[2,], along = 3)
    }
    
    pars <- lapply(sampler$samples, gibbs_step, sampler, prior, hyper)
    if(n_cores > 1){
      proposals=mclapply(X=1:pmwgs$n_subjects,FUN = new_particle, data, particles, pars, eff_mu, 
                         eff_var, mix, pmwgs$ll_func, epsilon, subjects, n_chains, mc.cores =n_cores)
    } else{
      proposals=lapply(X=1:pmwgs$n_subjects, FUN = new_particle, data, particles, pars, eff_mu, 
                       eff_var, mix, pmwgs$ll_func, epsilon, subjects, n_chains)
    }
    proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects, pmwgs$n_chains))

    j <- start_iter + i
    pmwgs$samples <- lapply(X = 1:n_chains, fill_samples, 
                            samples = pmwgs$samples,
                            group_level = pars,
                            proposals = proposals,
                            epsilon = epsilon,
                            j = j)
    

    
    # if(!is.null(pstar)){
    #   if(j > n0){
    #     acc <- alpha[1,] != pmwgs$samples$alpha[1,,(j-1)]
    #     epsilon<-max(update.epsilon(epsilon^2, acc, pstar, j, pmwgs$n_pars, alphaStar), epsilon_upper_bound)
    #   }
    # }
    # 
    # pmwgs$samples$epsilon[,j] <- epsilon
    
    if (stage == "adapt") {
      res <- test_sampler_adapted(pmwgs, n_unique, i, n_cores, get_conditionals)
      if (res == "success") {
        break
      } else if (res == "increase") {
        n_unique <- n_unique + .n_unique
      }
    }
  }
  if (display_progress) close(pb)
  if (stage == "adapt") {
    if (i == iter) {
      message(paste(
        "Particle Metropolis within Gibbs Sampler did not",
        "finish adaptation phase early (all", i, "iterations were",
        "run).\nYou should examine your samples and perhaps start",
        "a longer adaptation run."
      ))
    } else {
      pmwgs <- trim_na(pmwgs)
    }
  }
  return(pmwgs)
}

test_sampler_adapted <- function(pmwgs, n_unique, i, n_cores, conditionals_func) {
  n_pars <- length(pmwgs$par_names)
  if (i < n_unique) {
    return("continue")
  }
  test_samples <- extract_samples(pmwgs, stage = "adapt")
  # Only need to check uniqueness for one parameter
  first_par <- test_samples$alpha[1, , ]
  # Split the matrix into a list of vectors by subject
  # Needed for the case where every sample is unique for all subjects
  first_par_list <- split(first_par, seq(NROW(first_par)))
  # Get unique pars (new accepted particles) and check length for
  # all subjects is greater than unq_vals
  n_unique_sub <- lapply(lapply(first_par_list, unique), length)
  if (all(n_unique_sub > n_unique)) {
    message("Enough unique values detected: ", n_unique)
    message("Testing proposal distribution creation")
    attempt <- tryCatch({
      lapply(X = 1:pmwgs$n_subjects,FUN = conditionals_func,samples = test_samples, n_pars, i)
    },error=function(e) e, warning=function(w) w)
    if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
      warning("An problem was encountered creating proposal distribution")
      warning("Increasing required unique values and continuing adaptation")
      return("increase")
    }
    else {
      message("Successfully adapted after ", i, "iterations - stopping early")
      return("success")
    }
  }
  return("continue")
}


