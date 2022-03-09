### This is the pmwg sampler for standard puproses. It includes a multivariate normal with full covariance matrix on the group level. 

library(mvtnorm) ## For the multivariate normal.
library(MASS) ## For matrix inverse.
library(MCMCpack) #For inverse wishart
library(lme4)
library(parallel)
library(abind)
library(checkmate)
library(Rcpp)

source("multiGroup/utils_factors_nGroups_1Cov.R")
source("pmwg/messaging.R")

pmwgs <- function(data, pars, ll_func, n_factors, prior = NULL, constraintMat = NULL, par_idx = NULL) {
  ###Gets/sets priors, creates pmwgs object and stores essentials
  # Descriptives
  n_pars <- length(pars)
  subjects <- unique(data$subject)
  n_subjects <- length(subjects)
  groups <- unique(data$group)
  n_groups <- length(groups)
  parGroups <- unique(as.vector(par_idx))
  n_parGroups <- length(parGroups)
  group_idx <- aggregate(group ~ subject, data, mean)[,2]
  # Hyperparameters
  al <- 1 # hyperparameter shape for gamma on psi_inv
  bl <- 1/2 # hyperparameter rate for gamma on psi_inv
  nu <- 2 # hyperparameter shape for gamma on sig_err_inv 
  s2 <- 1/nu # hyperparameter rate for gamma on sig_err_inv
  
  # Storage for the samples.
  samples <- sample_store(pars, subjects, groups, n_factors)
  # Checking and default priors
  if (is.null(prior)) {
    lambda_idx <- c(c(1:n_factors), rep(n_factors, n_pars - n_factors))
    prior <- list(theta_mu_mean = rep(0, n_pars), 
                  theta_mu_var = rep(1, n_pars),
                  theta_lambda_var = 1)
  }
  if(is.null(par_idx)){
    par_idx <- matrix(0, nrow = n_pars, ncol = n_groups)
  }
  
  if(is.null(constraintMat)){
    constraintMat <- matrix(1, nrow = n_pars, ncol = n_factors)
    constraintMat[upper.tri(constraintMat)] <- 0 #Now you can't fix one of the diagonal values to 0
  }
  constraintMat <- constraintMat != 0 #For indexing
  if(any(diag(constraintMat) != 0)) signFix <- T
  
  
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- diag(1/prior$theta_mu_var) 
  prior$theta_lambda_invar <-1/prior$theta_lambda_var
  sampler <- list(
    data = data,
    par_names = pars,
    n_pars = n_pars,
    n_factors = n_factors,
    n_subjects = n_subjects,
    lambda_idx = lambda_idx,
    subjects = subjects,
    groups = groups,
    n_groups = n_groups,
    parGroups = parGroups,
    group_idx = group_idx,
    par_idx = par_idx,
    n_parGroups = n_parGroups,
    prior = prior,
    ll_func = ll_func,
    samples = samples,
    init = FALSE
  )
  #Hyper parameters
  attr(sampler, "al") <- al
  attr(sampler, "bl") <- bl
  attr(sampler, "nu") <- nu 
  attr(sampler, "s2") <- s2
  attr(sampler, "signFix") <- signFix
  attr(sampler, "constraintMat") <- constraintMat
  class(sampler) <- "pmwgs"
  sampler
}

init <- function(pmwgs, start_mu = NULL, start_var = NULL, start_psi_inv = NULL, 
                 start_sig_err_inv = NULL, start_lambda = NULL, 
                 display_progress = TRUE, particles = 1000, n_cores = 1, 
                 epsilon = NULL, useC = T) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  if (is.null(start_mu)) start_mu <- stats::rnorm(pmwgs$n_pars, sd = 1)
  #Even though we're using factor analysis, I'm using some start covariance matrix for startpoints to get reasonabe spread
  if (is.null(start_var)) {
    start_var <- MCMCpack::riwish(pmwgs$n_pars * 3,diag(pmwgs$n_pars))
  }
  
  # Factor analysis start points
  if(is.null(start_psi_inv)) start_psi_inv <- diag(1, pmwgs$n_factors)
  if(is.null(start_sig_err_inv)) start_sig_err_inv <- diag(1, pmwgs$n_pars)
  if(is.null(start_lambda)){
    start_lambda <- matrix(0, nrow = pmwgs$n_pars, ncol = pmwgs$n_factors)
    start_lambda[1:pmwgs$n_factors, 1:pmwgs$n_factors] <- diag(1, pmwgs$n_factors)
  }
  
  if(useC){
    sourceCpp("pmwg/utilityFunctions.cpp")
    dmv <<- dmvnrm_arma_fast
  } else{
    dmv <<- mvtnorm::dmvnorm
  }
  
  if(n_cores > 1){
    proposals <- mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals,start_mu = start_mu, 
                          start_var = start_var, n_particles = particles, pmwgs = pmwgs, mc.cores = n_cores)
  } else{
    proposals <- lapply(X=1:pmwgs$n_subjects,FUN=start_proposals,start_mu = start_mu, 
                        start_var = start_var, n_particles = particles, pmwgs = pmwgs)
  }
  proposals <- simplify2array(proposals)
  pmwgs$init <- TRUE
  
  pmwgs$samples$theta_lambda[,,1] <- start_lambda
  pmwgs$samples$theta_lambda_orig[,,1] <- start_lambda
  pmwgs$samples$theta_sig_err_inv[,,1] <- start_sig_err_inv
  pmwgs$samples$theta_psi_inv[,,1] <- start_psi_inv
  pmwgs$samples$theta_mu[,, 1] <- start_mu
  pmwgs$samples$theta_var[, , 1] <- start_var
  pmwgs$samples$theta_eta[,,1] <- matrix(0, nrow = pmwgs$n_subjects, ncol = pmwgs$n_factors)
  
  pmwgs$samples$alpha[, , 1] <- do.call(cbind, proposals[1,])
  pmwgs$samples$subj_ll[, 1] <- unlist(proposals[2,])
  pmwgs$samples$idx <- 1
  pmwgs$samples$epsilon[,1] <- rep(set_epsilon(pmwgs$n_pars, epsilon), pmwgs$n_subjects)
  pmwgs$samples$origin[,1] <- rep(2, pmwgs$n_subjects)
  return(pmwgs)
}

start_proposals <- function(s, start_mu, start_var, n_particles, pmwgs){
  #Draw the first start point
  proposals <- particle_draws(n_particles, start_mu, start_var)
  colnames(proposals) <- rownames(pmwgs$samples$theta_mu) # preserve par names
  lw <- apply(proposals,1,pmwgs$ll_func,data = pmwgs$data[pmwgs$data$subject == pmwgs$subjects[s], ])
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx]))
}

gibbs_step_factor <- function(sampler){
  # Gibbs step for group means with parameter expanded factor analysis from Ghosh & Dunson 2009
  # mu = theta_mu, var = theta_var
  last <- last_sample(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior
  
  #extract previous values (for ease of reading)
  
  alpha <- t(last$alpha)
  n_subjects <- sampler$n_subjects
  n_pars <- sampler$n_pars
  n_factors <- sampler$n_factors
  constraintMat <- hyper$constraintMat
  
  eta <- matrix(last$eta, n_subjects, n_factors)
  psi_inv <- matrix(last$psi_inv, n_factors)
  sig_err_inv <- last$sig_err_inv
  lambda <- matrix(last$lambda, n_pars, n_factors)
  mu <- last$mu
  alphatilde <- alpha
  #Update mu
  for(p in 1:sampler$n_parGroups){
    parGroup <- sampler$parGroups[p]
    par_idx <- unique(which(sampler$par_idx == parGroup, arr.ind = T)[,1])
    parGroup_num <- as.numeric(strsplit(parGroup,",")[[1]])
    group_idx <- which(sampler$group_idx %in% parGroup_num)
    mu_sig <- solve(length(group_idx) * sig_err_inv[par_idx, par_idx, drop = F] + prior$theta_mu_invar[par_idx, par_idx, drop = F])
    mu_mu <- mu_sig %*% (sig_err_inv[par_idx, par_idx, drop = F] %*% colSums(alpha[group_idx,par_idx, drop = F]
                         - eta[group_idx,,drop = F] %*% t(lambda[par_idx,, drop = F])) 
                         + prior$theta_mu_invar[par_idx, par_idx,drop = F] %*% prior$theta_mu_mean[par_idx])
    mu[par_idx,parGroup_num] <- rmvnorm(1, mu_mu, mu_sig)
  }
  for (g in 1:sampler$n_groups){
    group_idx <- which(sampler$group_idx == g)
    alphatilde[group_idx,] <- sweep(alphatilde[group_idx,], 2, mu[,g])
  }
  
  
  #Update eta, I do this one first since I don't want to save eta
  eta_sig <- solve(psi_inv + t(lambda) %*% sig_err_inv %*% lambda)
  eta_mu <- eta_sig %*% t(lambda) %*% sig_err_inv %*% t(alphatilde)
  eta[,] <- t(apply(eta_mu, 2, FUN = function(x){rmvnorm(1, x, eta_sig)}))
  
  #Update sig_err
  sig_err_inv <- diag(rgamma(n_pars,shape=(hyper$nu+n_subjects)/2, rate=(hyper$nu*hyper$s2+ colSums((alphatilde - eta %*% t(lambda))^2))/2))
  
  #Update lambda
  for (j in 1:n_pars) {
    constraint <- constraintMat[j,] #T if item is not constraint (bit confusing tbh)
    if(any(constraint)){ #Don't do this if there are no free entries in lambda
      etaS <- eta[,constraint]
      lambda_sig <- solve(sig_err_inv[j,j] * t(etaS) %*% etaS + prior$theta_lambda_invar * diag(1,sum(constraint)))
      lambda_mu <- (lambda_sig * sig_err_inv[j,j]) %*% (t(etaS) %*% alphatilde[,j])
      lambda[j,constraint] <- rmvnorm(1,lambda_mu,lambda_sig)
    }
  }
  
  #Update psi_inv
  psi_inv[,] <- diag(rgamma(n_factors ,shape=(hyper$al+n_subjects)/2,rate=hyper$bl+colSums(eta^2)/2), n_factors)
  
  lambda_orig <- lambda
  #If the diagonals of lambda aren't constrained to be 1, we should fix the signs
  if(hyper$signFix){
    for(l in 1:n_factors){
      mult <- ifelse(lambda[l, l] < 0, -1, 1) #definitely a more clever function for this 
      lambda_orig[,l] <- mult * lambda[, l]
    }
  }
  
  var <- lambda_orig %*% solve(psi_inv) %*% t(lambda_orig) + diag(1/diag((sig_err_inv)))
  return(list(tmu = mu, tvar = var, lambda = lambda, lambda_orig = lambda_orig, eta = eta,
              sig_err_inv = sig_err_inv, psi_inv = psi_inv, alpha = last$alpha))
}

new_particle <- function (s, data, num_particles, parameters, eff_mu = NULL, 
                          eff_var = NULL, mix_proportion = c(0.5, 0.5, 0), 
                          likelihood_func = NULL, epsilon = NULL, subjects, groups) 
{
  eff_mu <- eff_mu[, s]
  eff_var <- eff_var[, , s]
  mu <- parameters$tmu[,groups[s]]
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
  lp <- dmv(x = proposals, mean = mu, sigma = var, 
            log = TRUE)
  prop_density <- dmv(x = proposals, mean = subj_mu, 
                      sigma = var * (epsilon[s]^2))
  if (mix_proportion[3] == 0) {
    eff_density <- 0
  }
  else {
    eff_density <- dmv(x = proposals, mean = eff_mu, sigma = eff_var)
  }
  lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] * 
                                             prop_density) + (mix_proportion[3] * eff_density))
  l <- lw + lp - lm
  weights <- exp(l - max(l))
  idx <- sample(x = num_particles, size = 1, prob = weights)
  origin <- min(which(idx <= cumuNumbers))
  return(list(ll = lw[idx], origin = origin, proposal = proposals[idx, ]))
}

get_conditionals <- function(s, samples, n_pars, iteration, groups){
  pts2_unwound <- apply(samples$theta_var,3,unwind)
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu[,groups[s],],pts2_unwound)
  mu_tilde <- apply(all_samples, 1, mean)
  var_tilde <- stats::var(t(all_samples))
  condmvn <- condMVNorm::condMVN(mean = mu_tilde, sigma = var_tilde,
                                 dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                                 X.given = c(samples$theta_mu[,groups[s],iteration], unwind(samples$theta_var[,,iteration])))
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
                      pdist_update_n = ifelse(stage == "sample", 500, NA),
                      useC = T
) {
  # Set defaults for NULL values
  mix <- set_mix(stage, mix)
  # Set necessary local variables
  .n_unique <- n_unique
  # Set stable (fixed) new_sample argument for this run
  n_pars <- length(pmwgs$par_names)
  
  if(useC){
    sourceCpp("pmwg/utilityFunctions.cpp")
    dmv <<- dmvnrm_arma_fast
  } else{
    dmv <<- mvtnorm::dmvnorm
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
  
  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)
  # create progress bar
  if (display_progress) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx
  
  eff_mu <- NULL
  eff_var <- NULL
  data <- pmwgs$data
  subjects <- pmwgs$subjects
  group_idx <- pmwgs$group_idx
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
                              n_pars, iteration, group_idx, mc.cores = n_cores)
      } else{
        conditionals=lapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, 
                            n_pars, group_idx, iteration)
      }
      conditionals <- simplify2array(conditionals)
      eff_mu <- do.call(cbind, conditionals[1,])
      eff_var <- abind(conditionals[2,], along = 3)
    }
    
    pars <- gibbs_step_factor(pmwgs)
    if(n_cores > 1){
      proposals=mclapply(X=1:pmwgs$n_subjects,FUN = new_particle, data, particles, pars, eff_mu, 
                         eff_var, mix, pmwgs$ll_func, epsilon, subjects, group_idx, mc.cores =n_cores, mc.cleanup = T)
    } else{
      proposals=lapply(X=1:pmwgs$n_subjects, FUN = new_particle, data, particles, pars, eff_mu, 
                       eff_var, mix, pmwgs$ll_func, epsilon, subjects, group_idx)
    }
    proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects))
    ll <- proposals[1,]
    origin <- proposals[2,]
    alpha <- proposals[3:(pmwgs$n_pars+2),]
    
    j <- start_iter + i
    
    pmwgs$samples$theta_lambda[,,j] <- pars$lambda
    pmwgs$samples$theta_lambda_orig[,,j] <- pars$lambda_orig
    pmwgs$samples$theta_sig_err_inv[,,j] <- pars$sig_err_inv
    pmwgs$samples$theta_psi_inv[,,j] <- pars$psi_inv
    pmwgs$samples$theta_mu[,, j] <- pars$tmu
    pmwgs$samples$theta_var[, , j] <- pars$tvar
    pmwgs$samples$theta_eta[, , j] <- pars$eta
    
    pmwgs$samples$alpha[, , j] <- alpha
    pmwgs$samples$idx <- j
    pmwgs$samples$subj_ll[, j] <- ll
    pmwgs$samples$origin[,j] <- origin
    
    if(!is.null(pstar)){
      if(j > n0){
        acc <- alpha[1,] != pmwgs$samples$alpha[1,,(j-1)]
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
      if(n_cores > 1){
        mclapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, n_pars, i, pmwgs$group_idx, mc.cores = n_cores)
      } else{
        lapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, n_pars, i, pmwgs$group_idx)
      }
    },error=function(e) e, warning=function(w) w)
    if (any(class(attempt) %in% "warning")) {
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


