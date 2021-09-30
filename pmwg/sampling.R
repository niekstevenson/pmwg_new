library(mvtnorm) ## For the multivariate normal.
library(MASS) ## For matrix inverse.
library(MCMCpack) #For inverse wishart
library(lme4)
library(parallel)
library(abind)
library(checkmate)
library(Rcpp)

source("pmwg/utils.R")
source("pmwg/messaging.R")

init <- function(pmwgs, start_mu = NULL, start_sig = NULL,
         display_progress = TRUE, particles = 1000, n_cores = 1, epsilon = NULL, useC = T) {
  # If no starting point for group mean just use zeros
  if (is.null(start_mu)) start_mu <- stats::rnorm(pmwgs$n_pars, sd = 1)
  # If no starting point for group var just sample from inverse wishart
  if (is.null(start_sig)) {
    start_sig <- MCMCpack::riwish(pmwgs$n_pars * 3,diag(pmwgs$n_pars))
  }
  
  if(useC){
    sourceCpp("pmwg/utilityFunctions.cpp")
    rmv <<- mvrnorm_arma
    dmv <<- dmvnrm_arma_fast
  } else{
    rmv <<- mvtnorm::rmvnorm
    dmv <<- mvtnorm::dmvnorm
  }
  
  # Sample the mixture variables' initial values.
  a_half <- 1 / stats::rgamma(n = pmwgs$n_pars, shape = 0.5, scale = 1)
  # Create and fill initial random effects for each subject
  alpha <- array(NA, dim = c(pmwgs$n_pars, pmwgs$n_subjects))
  likelihoods <- array(NA_real_, dim = c(pmwgs$n_subjects))
  if(n_cores > 1){
    proposals <- mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals,start_mu = start_mu, 
                        start_sig = start_sig, n_particles = particles, pmwgs = pmwgs, mc.cores = n_cores)
  } else{
    proposals <- lapply(X=1:pmwgs$n_subjects,FUN=start_proposals,start_mu = start_mu, 
               start_sig = start_sig, n_particles = particles, pmwgs = pmwgs)
  }
  proposals <- simplify2array(proposals)
  pmwgs$init <- TRUE
  pmwgs$samples$theta_mu[, 1] <- start_mu
  pmwgs$samples$theta_sig[, , 1] <- start_sig
  pmwgs$samples$alpha[, , 1] <- do.call(cbind, proposals[1,])
  pmwgs$samples$last_theta_sig_inv <- MASS::ginv(start_sig)
  pmwgs$samples$subj_ll[, 1] <- unlist(proposals[2,])
  pmwgs$samples$a_half[, 1] <- a_half
  pmwgs$samples$idx <- 1
  pmwgs$samples$epsilon[,1] <- rep(set_epsilon(pmwgs$n_pars, epsilon), pmwgs$n_subjects)
  pmwgs$samples$origin[,1] <- rep(2, pmwgs$n_subjects)
  return(pmwgs)
}

start_proposals <- function(s, start_mu, start_sig, n_particles, pmwgs){
  proposals <- particle_draws(n_particles, start_mu, start_sig)
  colnames(proposals) <- rownames(pmwgs$samples$theta_mu) # preserve par names
  lw <- apply(proposals,1,pmwgs$ll_func,data = pmwgs$data[pmwgs$data$subject == pmwgs$subjects[s], ])
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx]))
}

gibbs_step <- function(sampler){
  # Get single iter versions, tmu = theta_mu, tsig = theta_sig
  last <- last_sample(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior
  
  # Here mu is group mean, so we are getting mean and variance
  var_mu <- MASS::ginv(sampler$n_subjects * last$tsinv + prior$theta_mu_invar)
  mean_mu <- as.vector(var_mu %*% (last$tsinv %*% apply(last$alpha, 1, sum) +
                                     prior$theta_mu_invar %*% prior$theta_mu_mean))
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  # New sample for mu.
  tmu <- mvtnorm::rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
  names(tmu) <- sampler$par_names
  
  # New values for group var
  theta_temp <- last$alpha - tmu
  cov_temp <- (theta_temp) %*% (t(theta_temp))
  B_half <- 2 * hyper$v_half * diag(1 / last$a_half) + cov_temp # nolint
  tsig <- MCMCpack::riwish(hyper$k_half, B_half) # New sample for group variance
  tsinv <- MASS::ginv(tsig)
  
  # Sample new mixing weights.
  a_half <- 1 / stats::rgamma(n = sampler$n_pars,shape = hyper$v_shape,
                              scale = 1 / (hyper$v_half + diag(tsinv) + hyper$A_half))
  return(list(tmu = tmu,tsig = tsig,tsinv = tsinv,a_half = a_half,alpha = last$alpha))
}

new_particle <- function (s, data, num_particles, parameters, eff_mu = NULL, 
                          eff_sig2 = NULL, mix_proportion = c(0.5, 0.5, 0), 
                          likelihood_func = NULL, epsilon = NULL) 
{
  eff_mu <- eff_mu[, s]
  eff_sig2 <- eff_sig2[, , s]
  mu <- parameters$tmu
  sig2 <- parameters$tsig
  subj_mu <- parameters$alpha[, s]
  particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
  cumuNumbers <- cumsum(particle_numbers)
  pop_particles <- particle_draws(particle_numbers[1], mu, 
                                  sig2)
  ind_particles <- particle_draws(particle_numbers[2], subj_mu, 
                                  sig2 * epsilon[s]^2)
  if(mix_proportion[3] == 0){
    eff_particles <- NULL
  } else{
    eff_particles <- particle_draws(particle_numbers[3], eff_mu, eff_sig2)
  }
  proposals <- rbind(pop_particles, ind_particles, eff_particles)
  colnames(proposals) <- names(mu)
  proposals[1, ] <- subj_mu
  lw <- apply(proposals, 1, likelihood_func, data = data[data$subject==s,])
  lp <- dmv(x = proposals, mean = mu, sigma = sig2, 
                         log = TRUE)
  prop_density <- dmv(x = proposals, mean = subj_mu, 
                                   sigma = sig2 * (epsilon[s]^2))
  if (mix_proportion[3] == 0) {
    eff_density <- 0
  }
  else {
    eff_density <- dmv(x = proposals, mean = eff_mu, sigma = eff_sig2)
  }
  lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] * 
                                             prop_density) + (mix_proportion[3] * eff_density))
  l <- lw + lp - lm
  weights <- exp(l - max(l))
  idx <- sample(x = num_particles, size = 1, prob = weights)
  origin <- min(which(idx <= cumuNumbers))
  return(list(proposal = proposals[idx, ], ll = lw[idx], origin = origin))
}
get_conditionals <- function(s, samples, n_pars, iteration){
  pts2_unwound <- apply(samples$theta_sig,3,unwind)
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu,pts2_unwound)
  mu_tilde <- apply(all_samples, 1, mean)
  sigma_tilde <- stats::var(t(all_samples))
  condmvn <- condMVNorm::condMVN(mean = mu_tilde, sigma = sigma_tilde,
                                 dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                                 X.given = c(samples$theta_mu[,iteration], unwind(samples$theta_sig[,,iteration])))
  return(list(eff_mu = condmvn$condMean, eff_sig2 = condmvn$condVar))
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
                      pdist_update_n = ifelse(stage == "sample", 500, NA)) {
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
  
  alphaStar=-qnorm(pstar/2) #Idk about this one
  n0=round(5/(pstar*(1-pstar))) #Also not questioning this math for now
  if(is.null(epsilon)){
    epsilon <- pmwgs$samples$epsilon[,ncol(pmwgs$samples$epsilon)]
  }
  
  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)
  # create progress bar
  if (display_progress) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx
  
  eff_mu <- NULL
  eff_sig2 <- NULL
  data <- pmwgs$data
  # Main iteration loop
  for (i in 1:iter) {
    accRate <- mean(accept_rate(pmwgs))
    if (display_progress) {
      update_progress_bar(pb, i, extra = accRate)
    }
    # Create/update efficient proposal distribution if we are in sampling phase.
    if(stage == "sample" & (i %% pdist_update_n == 0 || i == 1)){
      test_samples <- extract_samples(pmwgs, stage = c("adapt", "sample"))
      iteration <- dim(test_samples$theta_mu)[2]
      if(n_cores > 1){
        conditionals=mclapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, 
                              n_pars, iteration, mc.cores = n_cores)
      } else{
        conditionals=lapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, 
                            n_pars, iteration)
      }
      conditionals <- simplify2array(conditionals)
      eff_mu <- do.call(cbind, conditionals[1,])
      eff_sig2 <- abind(conditionals[2,], along = 3)
    }
    
    pars <- gibbs_step(pmwgs)
    if(n_cores > 1){
      proposals=mclapply(X=1:pmwgs$n_subjects,FUN = new_particle, data, particles, pars, eff_mu, 
                   eff_sig2, mix, pmwgs$ll_func, epsilon, mc.cores =n_cores, mc.cleanup = T)
    } else{
      proposals=lapply(X=1:pmwgs$n_subjects, FUN = new_particle, data, particles, pars, eff_mu, 
                 eff_sig2, mix, pmwgs$ll_func, epsilon)
    }
    proposals <- simplify2array(proposals)
    alpha <- do.call(cbind, proposals[1,])
    ll <- unlist(proposals[2,])
    origin <- unlist(proposals[3,])

    j <- start_iter + i
    
    pmwgs$samples$theta_mu[, j] <- pars$tmu
    pmwgs$samples$theta_sig[, , j] <- pars$tsig
    pmwgs$samples$last_theta_sig_inv <- pars$tsinv
    pmwgs$samples$alpha[, , j] <- alpha
    pmwgs$samples$idx <- j
    pmwgs$samples$subj_ll[, j] <- ll
    pmwgs$samples$origin[,j] <- origin
    pmwgs$samples$a_half[, j] <- pars$a_half
    
    if(!is.null(pstar)){
      if(j > n0){
        acc <- pmwgs$samples$alpha[1,,j] != pmwgs$samples$alpha[1,,(j-1)]
        epsilon<-update.epsilon(epsilon^2, acc, pstar, j, pmwgs$n_pars, alphaStar)
      }
    }
    
    pmwgs$samples$epsilon[,j] <- epsilon
    
    if (stage == "adapt") {
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

update.epsilon<- function(epsilon2, acc, p, i, d, alpha) {
  c=((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
  Theta=log(sqrt(epsilon2))
  Theta=Theta+c*(acc-p)/max(200, i/d)
  return(exp(Theta))
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
    attempt <- try({
      if(n_cores > 1){
        mclapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, n_pars, i, mc.cores = n_cores)
      } else{
        lapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, n_pars, i)
      }
    })
    if (class(attempt) == "try-error") {
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

particle_draws <- function(n, mu, covar) {
  if (n <= 0) {
    return(NULL)
  }
  rmv(n, mu, covar)
}
