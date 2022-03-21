### This is the pmwg sampler for standard puproses. It includes a multivariate normal with full covariance matrix on the group level. 

require(MASS) ## For matrix inverse.
require(MCMCpack) #For inverse wishart
require(lme4)
require(parallel)
require(mvtnorm) ## For the multivariate normal.
require(condMVNorm)
require(magic)
require(abind)

source("pmwg/messaging.R")

pmwgs <- function(data, pars, ll_func, prior = NULL, ...) {
  # Storage for the samples.
  subjects <- unique(data$subject)
  samples <- sample_store(pars, subjects, ...)
  sampler <- list(
    data = data,
    par_names = pars,
    subjects = subjects,
    n_pars = length(pars),
    n_subjects = length(subjects),
    ll_func = ll_func,
    samples = samples,
    init = FALSE
  )
  class(sampler) <- "pmwgs"
  sampler <- add_info(sampler, prior, ...)
  return(sampler)
}

init <- function(pmwgs, start_mu = NULL, start_var = NULL,
                 verbose = FALSE, particles = 1000, n_cores = 1, epsilon = NULL) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  startpoints <- get_startpoints(pmwgs, start_mu, start_var)
  if(n_cores > 1){
    proposals <- mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals,start_mu = startpoints$tmu, 
                          start_var = startpoints$tvar, n_particles = particles, pmwgs = pmwgs, mc.cores = n_cores)
  } else{
    proposals <- lapply(X=1:pmwgs$n_subjects,FUN=start_proposals,start_mu = startpoints$tmu, 
                        start_var = startpoints$tvar, n_particles = particles, pmwgs = pmwgs)
  }
  proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects))
  
  # Sample the mixture variables' initial values.
  if(is.null(epsilon)) epsilon <- rep(set_epsilon(pmwgs$n_pars, verbose), pmwgs$n_subjects)
  if(length(epsilon) == 1) epsilon <- rep(epsilon, pmwgs$n_subjects)
  pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = startpoints, proposals = proposals,
                                epsilon = epsilon, j = 1, n_pars = pmwgs$n_pars)
  pmwgs$init <- TRUE
  return(pmwgs)
}

fill_samples_base <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  samples$theta_mu[, j] <- group_level$tmu
  samples$theta_var[, , j] <- group_level$tvar
  samples$alpha[, , j] <- proposals[1:n_pars,]
  samples$subj_ll[, j] <- proposals[n_pars + 1,]
  samples$origin[,j] <- proposals[n_pars + 2,]
  samples$idx <- j
  samples$epsilon[,j] <- epsilon
  return(samples)
}


start_proposals <- function(s, start_mu, start_var, n_particles, pmwgs){
  #Draw the first start point
  proposals <- particle_draws(n_particles, start_mu, start_var)
  colnames(proposals) <- rownames(pmwgs$samples$theta_mu) # preserve par names
  lw <- apply(proposals,1,pmwgs$ll_func,data = pmwgs$data[pmwgs$data$subject == pmwgs$subjects[s], ])
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx], origin = 2))
}

new_particle <- function (s, data, num_particles, parameters, eff_mu = NULL, 
                          eff_var = NULL, mix_proportion = c(0.5, 0.5, 0), 
                          likelihood_func = NULL, epsilon = NULL, subjects) 
{
  eff_mu <- eff_mu[, s]
  eff_var <- eff_var[, , s]
  mu <- parameters$tmu
  var <- parameters$tvar
  subj_mu <- parameters$alpha[, s]
  particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
  cumuNumbers <- cumsum(particle_numbers)
  pop_particles <- particle_draws(particle_numbers[1], mu, 
                                  var)
  ind_particles <- particle_draws(particle_numbers[2], subj_mu, 
                                  var * epsilon[s]^2)
  if(mix_proportion[3] == 0){
    eff_particles <- NULL
  } else{
    eff_particles <- particle_draws(particle_numbers[3], eff_mu, eff_var)
  }
  proposals <- rbind(pop_particles, ind_particles, eff_particles)
  colnames(proposals) <- names(mu)
  proposals[1, ] <- subj_mu
  lw <- apply(proposals, 1, likelihood_func, data = data[data$subject==subjects[s],])
  lp <- dmvnorm(x = proposals, mean = mu, sigma = var, 
                log = TRUE)
  prop_density <- dmvnorm(x = proposals, mean = subj_mu, 
                          sigma = var * (epsilon[s]^2))
  if (mix_proportion[3] == 0) {
    eff_density <- 0
  }
  else {
    eff_density <- dmvnorm(x = proposals, mean = eff_mu, sigma = eff_var)
  }
  lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] * 
                                             prop_density) + (mix_proportion[3] * eff_density))
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
                      n_unique = ifelse(stage == "adapt", 20, NA),
                      min_unique = ifelse(stage == "adapt", 200, NA),
                      epsilon = NULL,
                      pstar = NULL,
                      mix = NULL,
                      pdist_update_n = ifelse(stage == "sample", 50, NA),
                      epsilon_upper_bound = 2,
                      n_cores_conditional = 1) {
  # Set defaults for NULL values
  mix <- set_mix(stage, mix, verbose)
  # Set necessary local variables
  .n_unique <- n_unique
  n_unique <- min_unique
  # Set stable (fixed) new_sample argument for this run
  n_pars <- length(pmwgs$par_names)
  
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
  
  eff_mu <- NULL
  eff_var <- NULL
  data <- pmwgs$data
  subjects <- pmwgs$subjects
  # Main iteration loop
  for (i in 1:iter) {
    accRate <- mean(accept_rate(pmwgs))
    if (verbose) {
      update_progress_bar(pb, i, extra = accRate)
    }
    # Create/update efficient proposal distribution if we are in sampling phase.
    if(stage == "sample" & (i %% pdist_update_n == 0 || i == 1)){
      test_samples <- extract_samples(pmwgs, stage = c("adapt", "sample"))
      iteration <- dim(test_samples$theta_mu)[2]
      if(n_cores_conditional > 1){
        conditionals=mclapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, 
                              n_pars, iteration, mc.cores = n_cores_conditional)
      } else{
        conditionals=lapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, 
                            n_pars, iteration)
      }
      conditionals <- simplify2array(conditionals)
      eff_mu <- do.call(cbind, conditionals[1,])
      eff_var <- abind(conditionals[2,], along = 3)
    }
    
    pars <- gibbs_step(pmwgs)
    if(n_cores > 1){
      proposals=mclapply(X=1:pmwgs$n_subjects,FUN = new_particle, data, particles, pars, eff_mu, 
                         eff_var, mix, pmwgs$ll_func, epsilon, subjects, mc.cores =n_cores)
    } else{
      proposals=lapply(X=1:pmwgs$n_subjects, FUN = new_particle, data, particles, pars, eff_mu, 
                       eff_var, mix, pmwgs$ll_func, epsilon, subjects)
    }
    proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects))
    
    
    
    #Fill samples
    j <- start_iter + i
    pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = pars,
                                  proposals = proposals, epsilon = epsilon, j = j, n_pars = pmwgs$n_pars)
    
    #Update epsilon
    if(!is.null(pstar)){
      if(j > n0){
        acc <-  pmwgs$samples$alpha[1,,j] != pmwgs$samples$alpha[1,,(j-1)]
        epsilon<-pmin(update.epsilon(epsilon^2, acc, pstar, j, pmwgs$n_pars, alphaStar), epsilon_upper_bound)
      }
    }
    
    if (stage == "adapt") {
      res <- test_sampler_adapted(pmwgs, n_unique, i, n_cores_conditional, get_conditionals, verbose)
      if (res == "success") {
        break
      } else if (res == "increase") {
        n_unique <- n_unique + .n_unique
      }
    }
  }
  if (verbose) close(pb)
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

test_sampler_adapted <- function(pmwgs, n_unique, i, n_cores_conditional, conditionals_func, verbose = T) {
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
    if(verbose){
      message("Enough unique values detected: ", n_unique)
      message("Testing proposal distribution creation")
    }
    attempt <- tryCatch({
      if(n_cores_conditional > 1){
        mclapply(X = 1:pmwgs$n_subjects,FUN = conditionals_func,samples = test_samples, 
                 n_pars, i, mc.cores = n_cores_conditional)
      } else{
        lapply(X = 1:pmwgs$n_subjects,FUN = conditionals_func,samples = test_samples, 
               n_pars, i)
      }
    },error=function(e) e, warning=function(w) w)
    if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
      if(verbose){
        warning("A problem was encountered creating proposal distribution")
        warning("Increasing required unique values and continuing adaptation")
      }
      return("increase")
    }
    else {
      if(verbose) message("Successfully adapted after ", i, "iterations - stopping early")
      return("success")
    }
  }
  return("continue")
}

particle_draws <- function(n, mu, covar) {
  if (n <= 0) {
    return(NULL)
  }
  return(rmvnorm(n, mu, covar))
}

update.epsilon<- function(epsilon2, acc, p, i, d, alpha) {
  c=((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
  Theta=log(sqrt(epsilon2))
  Theta=Theta+c*(acc-p)/max(200, i/d)
  return(exp(Theta))
}

set_epsilon <- function(n_pars, verbose = T) {
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

set_mix <- function(stage, mix, verbose) {
  if (stage == "sample") {
    mix <- c(0.1, 0.2, 0.7)
  } else {
    mix <- c(0.5, 0.5, 0.0)
  }
  if(verbose) message(sprintf("mix has been set to c(%s) based on the stage being run",  paste(mix, collapse = ", ")))
  return(mix)
}

numbers_from_proportion <- function(mix_proportion, num_particles = 1000) {
  numbers <- stats::rbinom(n = 2, size = num_particles, prob = mix_proportion)
  if (mix_proportion[3] == 0) {
    numbers[3] <- 0
    numbers[2] <- num_particles - numbers[1]
  } else {
    numbers[3] <- num_particles - sum(numbers)
  }
  return(numbers)
}

extract_samples <- function(sampler, stage = c("adapt", "sample")) {
  samples <- sampler$samples
  stage_filter <- samples$stage %in% stage
  sampled_filter <- seq_along(samples$stage) <= samples$idx
  
  list(
    theta_mu = samples$theta_mu[, stage_filter & sampled_filter],
    theta_var = samples$theta_var[, , stage_filter & sampled_filter],
    alpha = samples$alpha[, , stage_filter & sampled_filter]
  )
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

trim_na <- function(sampler) {
  idx <- sampler$samples$idx
  sampler$samples$stage <- sampler$samples$stage[1:idx]
  for(obj_name in names(sampler$samples)){
    obj <- sampler$samples[[obj_name]]
    dimensions <- length(dim(obj))
    if(dimensions > 1 & obj_name != "last_theta_var_inv"){ #Ok not the cleanest
      if(dimensions == 2){ #There must be a cleaner way to to do this
        obj <- obj[,1:idx]
      } else{
        obj <- obj[,,1:idx]
      }
      sampler$samples[[obj_name]] <- obj
    }
  }
  return(sampler)
}

