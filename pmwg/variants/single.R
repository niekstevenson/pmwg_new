library(mvtnorm) ## For the multivariate normal.
library(MASS) ## For matrix inverse.
library(MCMCpack) #For inverse wishart
library(lme4)
library(parallel)

source("pmwg/messaging.R")

pmwgs <- function(data, pars, ll_func, prior = NULL) {
  # Descriptives
  n_pars <- length(pars)
  subjects <- unique(data$subject)
  n_subjects <- length(subjects)
  
  # Storage for the samples.
  samples <- sample_store(pars, subjects)
  # Checking and default priors
  if (is.null(prior)) {
    prior <- list(theta_mu_mean = rep(0, n_pars), theta_mu_var = rep(1, n_pars))
  }

  prior$theta_mu_invar <- 1/prior$theta_mu_var #Inverse of the matrix
  sampler <- list(
    data = data,
    par_names = pars,
    n_pars = n_pars,
    n_subjects = n_subjects,
    subjects = subjects,
    prior = prior,
    ll_func = ll_func,
    samples = samples,
    init = FALSE
  )
  sampler
}

sample_store <- function(par_names, subject_ids, iters = 1, stage = "init", ...) {
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  list(
    epsilon = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    origin = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    alpha = array(NA_real_,dim = c(n_pars, n_subjects, iters),dimnames = list(par_names, subject_ids, NULL)),
    stage = array(stage, iters),
    subj_ll = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL))
  )
}

init <- function(pmwgs, start_mu = NULL, start_var = NULL,
                        verbose = FALSE, particles = 1000, n_cores = 1, epsilon = NULL) {
  # Create and fill initial random effects for each subject
  likelihoods <- array(NA_real_, dim = c(pmwgs$n_subjects))
  if(n_cores > 1){
    proposals <- mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals, pmwgs$prior$theta_mu_mean, pmwgs$prior$theta_mu_var, 
                          particles, pmwgs, mc.cores = n_cores)
  } else{
    proposals <- lapply(X=1:pmwgs$n_subjects,FUN=start_proposals, pmwgs$prior$theta_mu_mean, pmwgs$prior$theta_mu_var,
                        particles, pmwgs = pmwgs)
  }
  
  if(is.null(epsilon)) epsilon <- rep(set_epsilon(pmwgs$n_pars, verbose), pmwgs$n_subjects)
  if(length(epsilon) == 1) epsilon <- rep(epsilon, pmwgs$n_subjects)
  
  proposals <- simplify2array(proposals)
  pmwgs$init <- TRUE
  pmwgs$samples$alpha[, , 1] <- do.call(cbind, proposals[1,])
  pmwgs$samples$subj_ll[, 1] <- unlist(proposals[2,])
  pmwgs$samples$origin[,1] <- 2
  pmwgs$samples$idx <- 1
  pmwgs$samples$epsilon[,1] <- epsilon
  return(pmwgs)
}

start_proposals <- function(s, start_mu, start_var, n_particles, pmwgs){
  proposals <- particle_draws(n_particles, start_mu, start_var, pmwgs$n_pars)
  colnames(proposals) <- pmwgs$par_names # preserve par names
  lw <- apply(proposals,1,pmwgs$ll_func,data = pmwgs$data[pmwgs$data$subject == pmwgs$subjects[s], ])
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx]))
}

new_particle_single <- function (s, data, num_particles, prior_mu, prior_sd, alpha, mix_proportion = c(0.05, 0.95), 
                                 likelihood_func = NULL, epsilon = NULL, n_pars, subjects) 
{
  subj_mu <- alpha[,s,]
  particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
  cumuNumbers <- cumsum(particle_numbers)
  
  particle_numbers[2] <- num_particles - particle_numbers[1]
  pop_particles <- particle_draws(particle_numbers[1], prior_mu, 
                                  prior_sd, n_pars)
  ind_particles <- particle_draws(particle_numbers[2], subj_mu, 
                                  prior_sd * epsilon[s]^2, n_pars)
  proposals <- rbind(pop_particles, ind_particles)
  colnames(proposals) <- names(subj_mu)
  proposals[1, ] <- subj_mu
  lw <- apply(proposals, 1, likelihood_func, data = data[data$subject==subjects[s],])
  lp <- rowSums(dnorm(x = proposals, mean = prior_mu, sd = prior_sd, 
                      log = TRUE))
  prop_density <- rowSums(dnorm(x = proposals, mean = subj_mu, 
                                sd = prior_sd * (epsilon[s]^2)))
  
  lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] * prop_density))
  l <- lw + lp - lm
  weights <- exp(l - max(l))
  idx <- sample(x = num_particles, size = 1, prob = weights)
  origin <- min(which(idx <= cumuNumbers))
  return(list(proposal = proposals[idx, ], ll = lw[idx], origin = origin))
}

run_stage <- function(pmwgs,
                             stage,
                             iter = 1000,
                             particles = 1000,
                             verbose = TRUE,
                             n_cores = 1,
                             epsilon = NULL,
                             pstar = NULL,
                             mix = NULL) {
  # Set defaults for NULL values
  mix <- c(0.05, 0.95)
  # Set stable (fixed) new_sample argument for this run
  n_pars <- length(pmwgs$par_names)
  
  if(pmwgs$n_subjects == 1) verbose = F
  # Display stage to screen
  if(verbose){
    msgs <- list(
      burn = "Phase 1: Burn in\n",
      adapt = "Phase 2: Adaptation\n",
      sample = "Phase 3: Sampling\n"
    )
    cat(msgs[[stage]])
  }

  
  alphaStar=-qnorm(pstar/2) #Idk about this one
  n0=round(5/(pstar*(1-pstar))) #Also not questioning this math for now
  if(is.null(epsilon)){
    epsilon <- pmwgs$samples$epsilon[,ncol(pmwgs$samples$epsilon)]
  }  
  if(length(epsilon) == 1){
    epsilon <- rep(epsilon, pmwgs$n_subjects)
  }
  
  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)
  # create progress bar
  if (verbose) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx
  
  prior_mu <- pmwgs$prior$theta_mu_mean
  prior_var <- pmwgs$prior$theta_mu_var
  data <- pmwgs$data
  subjects <- pmwgs$subjects
  n_subjects <- pmwgs$n_subjects
  # Main iteration loop
  for (i in 1:iter) {
    if (verbose) {
      update_progress_bar(pb, i, extra = mean(accept_rate(pmwgs)))
    }
    j <- start_iter + i
    if(n_cores > 1){
      proposals=mclapply(X=1:pmwgs$n_subjects,FUN = new_particle_single, data, particles, prior_mu, prior_var,
                         pmwgs$samples$alpha[,,j-1, drop = F], mix, pmwgs$ll_func, epsilon, n_pars, subjects, mc.cores = n_cores)
    } else{
      proposals=lapply(X=1:pmwgs$n_subjects,FUN = new_particle_single, data, particles, prior_mu, prior_var,
                       pmwgs$samples$alpha[,,j-1, drop = F], mix, pmwgs$ll_func, epsilon, n_pars, subjects)
    }
    proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects))
    alpha <- matrix(proposals[1:n_pars,], nrow = n_pars, ncol = n_subjects)
    ll <- proposals[n_pars + 1,]
    origin <- proposals[n_pars + 2,]
    
    pmwgs$samples$alpha[, , j] <- alpha
    pmwgs$samples$idx <- j
    pmwgs$samples$subj_ll[, j] <- ll
    pmwgs$samples$origin[,j] <- origin
    
    if(!is.null(pstar)){
      if(j > n0){
        acc <- alpha[1,] != pmwgs$samples$alpha[1,,(j-1)]
        epsilon <- epsilon + (acc - pstar)/((i+n0)*pstar*(1-pstar))
      }
    }
    
    pmwgs$samples$epsilon[,j] <- epsilon
  }
  if (verbose) close(pb)
  return(pmwgs)
}

particle_draws <- function(n, mu, sd, n_pars) {
  if (n <= 0) {
    return(NULL)
  }
  matrix(rnorm(n*n_pars, mean = mu, sd =sd), ncol = n_pars, byrow = T)
}

set_epsilon <- function(n_pars, verbose) {
  if (n_pars > 15) {
    epsilon <- 0.1
  } else if (n_pars > 10) {
    epsilon <- 0.3
  } else {
    epsilon <- 0.5
  }
  if(verbose) message(sprintf("Epsilon has been set to %.1f based on number of parameters",epsilon))
  return(epsilon)
}

numbers_from_proportion <- function(mix_proportion, num_particles = 1000) {
  numbers <- stats::rbinom(n = 2, size = num_particles, prob = mix_proportion)
  numbers[2] <- num_particles - numbers[1]
  return(numbers)
}

extend_sampler <- function(sampler, n_samples, stage) {
  # This function takes the sampler and extends it along the intended number of
  # iterations, to ensure that we're not constantly increasing our sampled object
  # by 1. 
  old <- sampler$samples
  par_names <- sampler$par_names
  subject_ids <- sampler$subjects
  start <- old$idx + 1
  end <- old$idx + n_samples
  for(obj_name in names(sampler$samples)){
    obj <- sampler$samples[[obj_name]]
    dimensions <- length(dim(obj))
    if(dimensions > 1 & !obj_name %in% c("last_theta_var_inv")){ #Ok not the cleanest
      new_dim <- c(rep(0, (dimensions -1)), n_samples)
      extension <- array(NA_real_, dim = dim(obj) +  new_dim, dimnames = dimnames(obj))
      if(dimensions == 2){ #There must be a cleaner way to to do this
        extension[, -(start:end)] <- obj
      } else{
        extension[,, -(start:end)] <- obj
      }
      sampler$samples[[obj_name]] <- extension
    }
  }
  sampler$samples$stage <- c(old$stage, rep(stage, n_samples))
  sampler
}

