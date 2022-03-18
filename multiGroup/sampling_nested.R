library(matrixcalc)
library(MASS) ## For matrix inverse.
library(MCMCpack) #For inverse wishart
library(lme4)
library(parallel)
library(mvtnorm) ## For the multivariate normal.
library(condMVNorm)

source("multiGroup/utils_nested.R")
source("pmwg/messaging.R")

pmwgs <- function(data, pars, ll_func, prior = NULL) {
  ###Gets/sets priors, creates pmwgs object and stores essentials
  # Descriptives
  n_pars <- length(pars)
  subjects <- unique(data$subject)
  groups <- unique(data$group)
  n_groups <- length(groups)
  subjectgroups <- aggregate(group ~ subject, data, mean)[,2] #!!!Doesn't work with lists or tibbles, should fix this later
  n_subjects <- length(subjects)
  # Hyperparameters
  v_half <- 2 # hyperparameter on Σ prior (Half-t degrees of freedom)
  A_half <- 1 # hyperparameter on Σ prior (Half-t scale) #nolint
  
  # Create all the output arrays for the samplers
  samples <- sample_store(pars, subjects, groups)
  # Checking and default priors
  if (is.null(prior)) {
    prior <- list(theta_mu_mean = rep(0, n_pars), theta_mu_var = diag(rep(1, n_pars)))
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- ginv(prior$theta_mu_var) #Inverse of the matrix

  sampler <- list(
    data = data,
    par_names = pars,
    n_pars = n_pars,
    n_subjects = n_subjects,
    subjects = subjects,
    groups = groups,
    subjectgroups = subjectgroups,
    n_groups = n_groups,
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

init <- function(pmwgs, start_mu = NULL, start_var = NULL,
                 display_progress = TRUE, particles = 1000, n_cores = 1, epsilon = NULL) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  if (is.null(start_mu)){
    start_mu <- rnorm(pmwgs$n_pars, sd = 1)
    start_mu_group <- rnorm(pmwgs$n_pars, sd = 1)
  }
  # If no starting point for group var just sample some
  if (is.null(start_var)) {
    start_var <- riwish(pmwgs$n_pars * 3,diag(pmwgs$n_pars))
    start_var_group <- riwish(pmwgs$n_pars * 3,diag(pmwgs$n_pars))
  }
  
  # Sample the mixture variables' initial values.
  a_half <- 1 / rgamma(n = pmwgs$n_pars, shape = 2, rate = 1)
  # Create and fill initial random effects for each subject
  likelihoods <- array(NA_real_, dim = c(pmwgs$n_subjects))
  if(n_cores > 1){
    proposals <- mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals,start_mu = start_mu, 
                          start_var = start_var, n_particles = particles, pmwgs = pmwgs, mc.cores = n_cores)
  } else{
    proposals <- lapply(X=1:pmwgs$n_subjects,FUN=start_proposals,start_mu = start_mu, 
                        start_var = start_var, n_particles = particles, pmwgs = pmwgs)
  }
  proposals <- simplify2array(proposals)
  pmwgs$init <- TRUE
  
  #Population level mean of the mean and variance of the mean.
  pmwgs$samples$theta_mu[, 1] <- start_mu
  pmwgs$samples$theta_var[, , 1] <- start_var
  #Population level mean of the mean and variance of the mean.
  pmwgs$samples$group_mu[,,1] <- start_mu
  pmwgs$samples$group_var[,,,1] <- start_var_group
  #Random effects
  pmwgs$samples$alpha[, , 1] <- do.call(cbind, proposals[1,])
  #Inverse of the variance of the means, use these in Gibbs steps. 
  pmwgs$samples$last_theta_var_inv <- ginv(start_var)
  pmwgs$samples$last_group_var_inv <-array(rep(ginv(start_var), pmwgs$n_groups), 
                                                 dim = c(pmwgs$n_pars, pmwgs$n_pars, pmwgs$n_groups))
  #Individual likelihoods
  pmwgs$samples$subj_ll[, 1] <- unlist(proposals[2,])
  #'mixing' weights for huang and wand priors (2013)
  pmwgs$samples$group_a_half[,,1] <- a_half
  pmwgs$samples$a_half[, 1] <- a_half
  pmwgs$samples$idx <- 1
  #Epsilon = scaling of the covariance matrix, origin is the location where the particle was accepted from (group, individual, or conditional 'efficient' mean)
  pmwgs$samples$epsilon[,1] <- rep(set_epsilon(pmwgs$n_pars, epsilon), pmwgs$n_subjects)
  pmwgs$samples$origin[,1] <- rep(2, pmwgs$n_subjects)
  return(pmwgs)
}

start_proposals <- function(s, start_mu, start_var, n_particles, pmwgs){
  #Draw the first start point
  proposals <- particle_draws(n_particles, start_mu, start_var)
  colnames(proposals) <- rownames(pmwgs$samples$theta_mu) # preserve par names
  #calculate their likelihoods
  lw <- apply(proposals,1,pmwgs$ll_func,data = pmwgs$data[pmwgs$data$subject == pmwgs$subjects[s], ])
  #Importance sampling
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx]))
}

gibbs_step_pop <- function(sampler){
  # Gibbs step for population means, with full covariance matrix estimation
  # tmu = theta_mu, tvar = theta_var
  last <- last_sample_nested(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior
  
  # Here mu is pop mean, so we are getting mean and variance of the mean
  var_mu <- ginv(sampler$n_groups * last$tvinv + prior$theta_mu_invar)
  mean_mu <- as.vector(var_mu %*% (last$tvinv %*% apply(last$group_mu, 1, sum) +
                                     prior$theta_mu_invar %*% prior$theta_mu_mean))
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  # New sample for mu.
  tmu <- rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
  names(tmu) <- sampler$par_names
  
  # New values for pop var (Huang and Wand 2013)
  theta_temp <- last$group_mu - tmu
  cov_temp <- (theta_temp) %*% (t(theta_temp))
  B_half <- 2 * hyper$v_half * diag(1 / last$a_half) + cov_temp # nolint
  tvar <- riwish(hyper$v_half + sampler$n_pars - 1 + sampler$n_groups, B_half) # New sample for group variance
  tvinv <- ginv(tvar)
  
  # Sample new mixing weights.
  a_half <- 1 / rgamma(n = sampler$n_pars,shape = (hyper$v_half + sampler$n_pars) / 2,
                              rate = hyper$v_half * diag(tvinv) + hyper$A_half)
  return(list(tmu = tmu,tvar = tvar,tvinv = tvinv,a_half = a_half,group_mu = last$group_mu))
}


gibbs_step_group <- function(group, sampler){
  # Gibbs step for group means, with full covariance matrix estimation
  # tmu = theta_mu, tvar = theta_var
  last <- last_sample_group(sampler$samples, group)
  hyper <- attributes(sampler)
  idx <- sampler$subjectgroups == group
  n_subjects <- sum(idx)
  
  # Here mu is group mean, so we are getting mean and variance
  var_mu <- ginv(n_subjects * last$tvinv + last$theta_mu_invar)
  mean_mu <- as.vector(var_mu %*% (last$tvinv %*% apply(last$alpha[,idx], 1, sum) +
                                     last$theta_mu_invar %*% last$theta_mu_mean))

  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  # New sample for mu.
  tmu <- rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
  names(tmu) <- sampler$par_names
  
  # New values for group var (Huang and Wand 2013)
  theta_temp <- last$alpha[,idx] - tmu
  cov_temp <- (theta_temp) %*% (t(theta_temp))
  B_half <- 2 * hyper$v_half * diag(1 / last$a_half) + cov_temp # nolint
  tvar <- riwish(hyper$v_half + sampler$n_pars - 1 + n_subjects, B_half) # New sample for group variance
  tvinv <- ginv(tvar)
  
  # Sample new mixing weights.
  a_half <- 1 / rgamma(n = sampler$n_pars, shape = (hyper$v_half + sampler$n_pars) / 2,
                              rate = hyper$v_half * diag(tvinv) + hyper$A_half)
  return(list(tmu = tmu,tvar = tvar,tvinv = tvinv,a_half = a_half, alpha = last$alpha))
}


new_particle <- function (s, data, num_particles, group_parameters, eff_mu = NULL, 
                          eff_var = NULL, mix_proportion = c(0.5, 0.5, 0), 
                          likelihood_func = NULL, epsilon = NULL, subjects, subjectgroups){
  #Select the group pars for the group this individual belonged to. 
  group_parameters <- group_parameters[[subjectgroups[s]]]
  #These are the efficient multivariate normals based on their individual particles
  eff_mu <- eff_mu[, s]
  eff_var <- eff_var[, , s]
  #Get mean and covariance matrix for that group
  group_mu <- group_parameters$tmu
  group_var <- group_parameters$tvar
  #Subject specific mean is based on their previous accepted particle
  subj_mu <- group_parameters$alpha[, s]
  #Get the particles division
  particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
  cumuNumbers <- cumsum(particle_numbers)
  #Sample from group and subject mean and if in 'sampling' stage also from efficcient mean and covariance
  group_particles <- particle_draws(particle_numbers[1], group_mu, 
                                    group_var)
  ind_particles <- particle_draws(particle_numbers[2], subj_mu, 
                                  group_var * epsilon[s]^2)
  if(mix_proportion[3] == 0){
    eff_particles <- NULL
  } else{
    eff_particles <- particle_draws(particle_numbers[3], eff_mu, eff_var)
  }
  #Group them all together
  proposals <- rbind(group_particles, ind_particles, eff_particles)
  colnames(proposals) <- names(group_mu)
  proposals[1, ] <- subj_mu #This still hurts rather than rbind
  #Do some maths on the particles likelihoods and their origins (are in the new estimation approaches pmwg paper)
  lw <- apply(proposals, 1, likelihood_func, data = data[data$subject==subjects[s],])
  lg <- dmvnorm(x = proposals, mean = group_mu, sigma = group_var, log = TRUE)
  prop_density <- dmvnorm(x = proposals, mean = subj_mu, sigma = group_var * (epsilon[s]^2))
  if (mix_proportion[3] == 0) {
    eff_density <- 0
  }
  else {
    eff_density <- dmvnorm(x = proposals, mean = eff_mu, sigma = eff_var)
  }
  #Importance sampling
  lm <- log(mix_proportion[1] * exp(lg) + mix_proportion[2] * prop_density + (mix_proportion[3] * eff_density))
  l <- lw + lg - lm
  weights <- exp(l - max(l))
  idx <- sample(x = num_particles, size = 1, prob = weights)
  origin <- min(which(idx <= cumuNumbers))
  return(list(proposal = proposals[idx, ], ll = lw[idx], origin = origin))
}

get_conditionals <- function(s, samples, n_pars, iteration, subjectgroups){
  # Create a conditional multivariate normal based on group mean and variance and individual mean
  group <- subjectgroups[s]
  group_var <- samples$group_var[,,group,]
  group_mu <- samples$group_mu[,group,]
  pts2_unwound <- apply(group_var,3,unwind)
  all_samples <- rbind(samples$alpha[, s,],group_mu,pts2_unwound)
  mu_tilde <- apply(all_samples, 1, mean)
  var_tilde <- var(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                                 dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                                 X.given = c(group_mu[,iteration], unwind(group_var[,,iteration])))
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
                      diag_only = F,
                      pdist_update_n = ifelse(stage == "sample", 500, NA)
) {
  # Set defaults for NULL values
  mix <- set_mix(stage, mix)
  # Set necessary local variables
  .n_unique <- n_unique
  # Set stable (fixed) new_sample argument for this run
  n_pars <- length(pmwgs$par_names)
  
  # Display stage to screen
  msgs <- list(
    burn = "Phase 1: Burn in\n",
    adapt = "Phase 2: Adaptation\n",
    sample = "Phase 3: Sampling\n"
  )
  cat(msgs[[stage]])
  
  alphaStar=-qnorm(pstar/2) #Idk about this one but necessary for algorithm proposed in Garthwaite, Yan and Sisson 2016 
  n0=round(5/(pstar*(1-pstar))) #Same here
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
  
  #Some constants
  eff_mu <- NULL
  eff_var <- NULL
  data <- pmwgs$data
  subjects <- pmwgs$subjects
  subjectgroups <- pmwgs$subjectgroups
  n_cores_group <- min(pmwgs$n_groups, n_cores)
  # Main iteration loop
  for (i in 1:iter) {
    accRate <- mean(accept_rate(pmwgs))
    if (display_progress) {
      update_progress_bar(pb, i, extra = accRate)
    }
    # Create/update efficient proposal distribution if we are in sampling phase.
    if(stage == "sample" & (i %% pdist_update_n == 0 || i == 1)){
      test_samples <- extract_samples(pmwgs, stage = c("adapt", "sample"))
      iteration <- dim(test_samples$group_mu)[3]
      if(n_cores > 1){
        conditionals=mclapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, 
                              n_pars, iteration, subjectgroups, mc.cores = n_cores)
      } else{
        conditionals=lapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, 
                            n_pars, iteration, subjectgroups)
      }
      conditionals <- simplify2array(conditionals)
      eff_mu <- do.call(cbind, conditionals[1,])
      eff_var <- abind(conditionals[2,], along = 3)
    }
    
    #Get parameters for pop mean through gibbs step
    pars_pop <- gibbs_step_pop(pmwgs)
    #Same but for group means
    if(n_cores_group > 1){
      pars_group=mclapply(X=1:pmwgs$n_groups,FUN = gibbs_step_group, pmwgs, mc.cores =n_cores_group, mc.cleanup = T)
    } else{
      pars_group=lapply(X=1:pmwgs$n_groups,FUN = gibbs_step_group, pmwgs)
    }
    #Do particle metropolis for all subjects
    if(n_cores > 1){
      proposals=mclapply(X=1:pmwgs$n_subjects,FUN = new_particle, data, particles, pars_group, eff_mu, 
                         eff_var, mix, pmwgs$ll_func, epsilon, subjects, subjectgroups, mc.cores =n_cores, mc.cleanup = T)
    } else{
      proposals=lapply(X=1:pmwgs$n_subjects, FUN = new_particle, data, particles, pars_group, eff_mu, 
                       eff_var, mix, pmwgs$ll_func, epsilon, subjects, subjectgroups)
    }
    
    proposals <- simplify2array(proposals)
    alpha <- do.call(cbind, proposals[1,])
    ll <- unlist(proposals[2,])
    origin <- unlist(proposals[3,])
    
    pars_group <- simplify2array(pars_group)
    
    j <- start_iter + i
    
    #Store all my output
    pmwgs$samples$group_mu[,, j] <- do.call(cbind, pars_group[1,])
    pmwgs$samples$group_var[,,,j] <- abind(pars_group[2,], along = 3)
    pmwgs$samples$last_group_var_inv <- abind(pars_group[3,], along = 3)
    pmwgs$samples$group_a_half[,, j] <- do.call(cbind, pars_group[4,])
    
    
    pmwgs$samples$theta_mu[, j] <- pars_pop$tmu
    pmwgs$samples$theta_var[, , j] <- pars_pop$tvar
    pmwgs$samples$last_theta_var_inv <- pars_pop$tvinv
    pmwgs$samples$alpha[, , j] <- alpha
    pmwgs$samples$idx <- j
    pmwgs$samples$subj_ll[, j] <- ll
    pmwgs$samples$origin[,j] <- origin
    pmwgs$samples$a_half[, j] <- pars_pop$a_half
    
    #Adaptive scaling of epsilon based on Garthwaite, Yan and Sisson 2016 
    if(!is.null(pstar)){
      if(j > n0){
        acc <- alpha[1,] != pmwgs$samples$alpha[1,,(j-1)]
        epsilon<-update.epsilon(epsilon^2, acc, pstar, j, pmwgs$n_pars, alphaStar)
      }
    }
    
    pmwgs$samples$epsilon[,j] <- epsilon
    
    if (stage == "adapt") {
      #If we're in adaptation stage check whether we can create a conditional
      res <- test_sampler_adapted(pmwgs, n_unique, i, n_cores)
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

test_sampler_adapted <- function(pmwgs, n_unique, i, n_cores) {
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
      lapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, n_pars, i,
               pmwgs$subjectgroups)
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


