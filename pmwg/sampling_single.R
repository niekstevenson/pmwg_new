library(mvtnorm) ## For the multivariate normal.
library(MASS) ## For matrix inverse.
library(MCMCpack) #For inverse wishart
library(lme4)
library(parallel)

source("pmwg/utils.R")
source("pmwg/messaging.R")

pmwgs <- function(data, pars, ll_func, prior = NULL) {
  ###Gets/sets priors, creates pmwgs object and stores essentials
  # Descriptives
  n_pars <- length(pars)
  subjects <- unique(data$subject)
  n_subjects <- length(subjects)
  
  # Hyperparameters
  v_half <- 2 # hyperparameter on Σ prior (Half-t degrees of freedom)
  A_half <- 1 # hyperparameter on Σ prior (Half-t scale) #nolint
  
  # Storage for the samples.
  samples <- sample_store(pars, subjects)
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
    init = FALSE
  )
  #Hyper parameters
  attr(sampler, "v_half") <- v_half
  attr(sampler, "A_half") <- A_half
  class(sampler) <- "pmwgs"
  sampler
}


init_single <- function(pmwgs, start_mu = NULL, start_var = NULL,
                        display_progress = TRUE, particles = 1000, n_cores = 1, epsilon = NULL) {
  # Create and fill initial random effects for each subject
  likelihoods <- array(NA_real_, dim = c(pmwgs$n_subjects))
  if(n_cores > 1){
    proposals <- mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals, pmwgs$prior$theta_mu_mean, pmwgs$prior$theta_mu_var, 
                          particles, pmwgs, mc.cores = n_cores)
  } else{
    proposals <- lapply(X=1:pmwgs$n_subjects,FUN=start_proposals, pmwgs$prior$theta_mu_mean, pmwgs$prior$theta_mu_var,
                        particles, pmwgs = pmwgs)
  }
  proposals <- simplify2array(proposals)
  pmwgs$init <- TRUE
  pmwgs$samples$alpha[, , 1] <- do.call(cbind, proposals[1,])
  pmwgs$samples$subj_ll[, 1] <- unlist(proposals[2,])
  pmwgs$samples$idx <- 1
  pmwgs$samples$epsilon[,1] <- rep(set_epsilon(pmwgs$n_pars, epsilon), pmwgs$n_subjects)
  return(pmwgs)
}

start_proposals <- function(s, start_mu, start_var, n_particles, pmwgs){
  proposals <- particle_draws(n_particles, start_mu, start_var)
  colnames(proposals) <- rownames(pmwgs$samples$theta_mu) # preserve par names
  lw <- apply(proposals,1,pmwgs$ll_func,data = pmwgs$data[pmwgs$data$subject == pmwgs$subjects[s], ])
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx]))
}

new_particle_single <- function (s, data, num_particles, prior_mu, prior_var, alpha, mix_proportion = c(0.05, 0.95), 
                                 likelihood_func = NULL, epsilon = NULL, subjects) 
{
  subj_mu <- alpha[, s]
  particle_numbers <- rbinom(n = 2, size = num_particles, prob = mix_proportion)
  particle_numbers[2] <- num_particles - particle_numbers[1]
  cumuNumbers <- cumsum(particle_numbers)
  pop_particles <- particle_draws(particle_numbers[1], prior_mu, 
                                prior_var)
  ind_particles <- particle_draws(particle_numbers[2], subj_mu, 
                                prior_var * epsilon[s]^2)
  proposals <- rbind(pop_particles, ind_particles)
  colnames(proposals) <- names(subj_mu)
  proposals[1, ] <- subj_mu
  lw <- apply(proposals, 1, likelihood_func, data = data[data$subject==subjects[s],])
  lp <- dmvnorm(x = proposals, mean = prior_mu, sigma = prior_var, 
                         log = TRUE)
  prop_density <- dmvnorm(x = proposals, mean = subj_mu, 
                                   sigma = prior_var * (epsilon[s]^2))
  
  lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] * prop_density))
  l <- lw + lp - lm
  weights <- exp(l - max(l))
  idx <- sample(x = num_particles, size = 1, prob = weights)
  origin <- min(which(idx <= cumuNumbers))
  return(list(proposal = proposals[idx, ], ll = lw[idx], origin = origin))
}

run_stage_single <- function(pmwgs,
                             stage,
                             iter = 1000,
                             particles = 1000,
                             display_progress = TRUE,
                             n_cores = 1,
                             epsilon = NULL,
                             pstar = NULL,
                             mix = NULL) {
  # Set defaults for NULL values
  mix <- c(0.05, 0.95)
  # Set stable (fixed) new_sample argument for this run
  n_pars <- length(pmwgs$par_names)
  
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
  
  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)
  # create progress bar
  if (display_progress) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx
  
  prior_mu <- pmwgs$prior$theta_mu_mean
  prior_var <- pmwgs$prior$theta_mu_var
  data <- pmwgs$data
  subjects <- pmwgs$subjects
  # Main iteration loop
  for (i in 1:iter) {
    if (display_progress) {
      update_progress_bar(pb, i, extra = mean(accept_rate(pmwgs)))
    }
    j <- start_iter + i
    if(n_cores > 1){
      proposals=mclapply(X=1:pmwgs$n_subjects,FUN = new_particle_single, data, particles, prior_mu, prior_var,
                         pmwgs$samples$alpha[,,j-1], mix, pmwgs$ll_func, epsilon, subjects, mc.cores = n_cores)
    } else{
      proposals=lapply(X=1:pmwgs$n_subjects,FUN = new_particle_single, data, particles, prior_mu, prior_var,
                       pmwgs$samples$alpha[,,j-1], mix, pmwgs$ll_func, epsilon, subjects)
    }
    proposals <- simplify2array(proposals)
    alpha <- do.call(cbind, proposals[1,])
    ll <- unlist(proposals[2,])
    origin <- unlist(proposals[3,])
    
    pmwgs$samples$alpha[, , j] <- alpha
    pmwgs$samples$idx <- j
    pmwgs$samples$subj_ll[, j] <- ll
    pmwgs$samples$origin[,j] <- origin
    
    if(!is.null(pstar)){
      if(j > n0){
        acc <- alpha[1,] != pmwgs$samples$alpha[1,,(j-1)]
        epsilon<-update.epsilon(epsilon^2, acc, pstar, j, n_pars, alphaStar)
      }
    }
    
    pmwgs$samples$epsilon[,j] <- epsilon
  }
  if (display_progress) close(pb)
  return(pmwgs)
}

update.epsilon<- function(epsilon2, acc, p, i, d, alpha) {
  c=((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
  Theta=log(sqrt(epsilon2))
  Theta=Theta+c*(acc-p)/max(200, i/d)
  return(exp(Theta))
}

particle_draws <- function(n, mu, covar) {
  if (n <= 0) {
    return(NULL)
  }
  rmvnorm(n, mu, covar)
}


