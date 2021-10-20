extract_samples <- function(sampler, stage = c("adapt", "sample")) {
  samples <- sampler$samples
  stage_filter <- samples$stage %in% stage
  sampled_filter <- seq_along(samples$stage) <= samples$idx
  
  list(
    theta_mu = samples$theta_mu[, stage_filter & sampled_filter],
    theta_sig = samples$theta_sig[, , stage_filter & sampled_filter],
    alpha = samples$alpha[, , stage_filter & sampled_filter]
  )
}

particle_draws <- function(n, mu, covar) {
  if (n <= 0) {
    return(NULL)
  }
  rmv(n, mu, covar)
}

update.epsilon<- function(epsilon2, acc, p, i, d, alpha) {
  c=((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
  Theta=log(sqrt(epsilon2))
  Theta=Theta+c*(acc-p)/max(200, i/d)
  return(exp(Theta))
}

unwind <- function(var_matrix, ...) {
  y <- t(chol(var_matrix))
  diag(y) <- log(diag(y))
  y[lower.tri(y, diag = TRUE)]
}

wind <- function(var_vector, ...) {
  n <- sqrt(2 * length(var_vector) + 0.25) - 0.5 ## Dim of matrix.
  if ((n * n + n) != (2 * length(var_vector))) stop("Wrong sizes in unwind.")
  out <- array(0, dim = c(n, n))
  out[lower.tri(out, diag = TRUE)] <- var_vector
  diag(out) <- exp(diag(out))
  out %*% t(out)
}


extend_sampler <- function(sampler, n_samples, stage) {
  old <- sampler$samples
  par_names <- sampler$par_names
  subject_ids <- sampler$subjects
  groups <- sampler$groups
  start <- old$idx + 1
  end <- old$idx + n_samples
  
  new_group_mu <- array(NA_real_,dim = dim(old$group_mu) + c(0, 0, n_samples), dimnames = list(par_names, groups, NULL))
  new_group_mu[,, - (start:end)] <- old$group_mu
  sampler$samples$group_mu <- new_group_mu
  
  new_group_sig <- array(NA_real_,dim = dim(old$group_sig) + c(0, 0, 0, n_samples),dimnames = list(par_names, par_names, groups, NULL))
  new_group_sig[,,, - (start:end)] <- old$group_sig
  sampler$samples$group_sig <- new_group_sig
  
  new_group_a_half <- array(NA_real_,dim = dim(old$group_a_half) + c(0,0, n_samples), dimnames = list(par_names, groups, NULL))
  new_group_a_half[,, - (start:end)] <- old$group_a_half
  sampler$samples$group_a_half <- new_group_a_half
  
  new_epsilon <- array(NA_real_,dim = dim(old$epsilon) + c(0, n_samples),dimnames = list(subject_ids, NULL))
  new_epsilon[, - (start:end)] <- old$epsilon
  sampler$samples$epsilon <- new_epsilon
  
  new_orig <- array(NA_real_,dim = dim(old$origin) + c(0, n_samples),dimnames = list(subject_ids, NULL))
  new_orig[, - (start:end)] <- old$origin
  sampler$samples$origin <- new_orig
  
  new_tmu <- array(NA_real_,dim = dim(old$theta_mu) + c(0, n_samples), dimnames = list(par_names, NULL))
  new_tmu[, - (start:end)] <- old$theta_mu
  sampler$samples$theta_mu <- new_tmu
  
  new_tsig <- array(NA_real_,dim = dim(old$theta_sig) + c(0, 0, n_samples),dimnames = list(par_names, par_names, NULL))
  new_tsig[, , - (start:end)] <- old$theta_sig
  sampler$samples$theta_sig <- new_tsig
  
  new_alph <- array(NA_real_,dim = dim(old$alpha) + c(0, 0, n_samples),dimnames = list(par_names, subject_ids, NULL))
  new_alph[, , - (start:end)] <- old$alpha
  sampler$samples$alpha <- new_alph
  
  new_sll <- array(NA_real_,dim = dim(old$subj_ll) + c(0, n_samples),dimnames = list(subject_ids, NULL))
  new_sll[, - (start:end)] <- old$subj_ll
  sampler$samples$subj_ll <- new_sll
  
  new_ahalf <- array(NA_real_,dim = dim(old$a_half) + c(0, n_samples), dimnames = list(par_names, NULL))
  new_ahalf[, - (start:end)] <- old$a_half
  sampler$samples$a_half <- new_ahalf
  
  sampler$samples$stage <- c(old$stage, rep(stage, n_samples))
  sampler
}

trim_na <- function(sampler) {
  idx <- sampler$samples$idx
  sampler$samples$theta_mu <- sampler$samples$theta_mu[, 1:idx]
  sampler$samples$theta_sig <- sampler$samples$theta_sig[, , 1:idx]
  sampler$samples$alpha <- sampler$samples$alpha[, , 1:idx]
  sampler$samples$subj_ll <- sampler$samples$subj_ll[, 1:idx]
  sampler$samples$a_half <- sampler$samples$a_half[, 1:idx]
  sampler$samples$stage <- sampler$samples$stage[1:idx]
  sampler$samples$epsilon <- sampler$samples$epsilon[,1:idx]
  sampler$samples$origin <- sampler$samples$origin[, 1:idx]
  sampler
}
relabel_samples <- function(sampler, indices, from="burn", to="adapt") {
  old_stage <- sampler$samples$stage
  if (!all(old_stage[indices] %in% from)) {
    stop(paste("Not all samples were from the", from, "stage"))
  }
  sampler$samples$stage[indices] <- to
  sampler
}


as_mcmc <- function(sampler, selection = "theta_mu", filter = stages) {
  stages <- c("burn", "adapt", "sample")
  if (all(filter %in% stages)) {
    filter <- which(sampler$samples$stage %in% filter)
  } else if (!all(filter %in% 1:sampler$samples$idx)) {
    stop("filter is not a vector of stage names, or integer vector of indices")
  }
  
  if (selection == "theta_mu") {
    return(coda::mcmc(t(sampler$samples$theta_mu[, filter])))
  } else if (selection == "theta_sig") {
    tsig <- sampler$samples$theta_sig[, , filter]
    return(stats::setNames(lapply(
      seq(dim(tsig)[1]),
      function(x) {
        coda::mcmc(t(tsig[x, , ]))
      }
    ), sampler$par_names))
  } else if (selection == "alpha") {
    alpha <- sampler$samples$alpha[, , filter]
    return(stats::setNames(lapply(
      seq(dim(alpha)[2]),
      function(x) {
        coda::mcmc(t(alpha[, x, ]))
      }
    ), sampler$subjects))
  }
  stop("Argument `selection` should be one of theta_mu, theta_sig, alpha")
}

sample_store <- function(par_names, subject_ids, groups, iters = 1, stage = "init") {
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  n_groups <- length(groups)
  list(
    epsilon = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    origin = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    alpha = array(NA_real_,dim = c(n_pars, n_subjects, iters),dimnames = list(par_names, subject_ids, NULL)),
    group_mu = array(NA_real_,dim = c(n_pars, n_groups, iters), dimnames = list(par_names, groups, NULL)),
    group_sig = array(NA_real_,dim = c(n_pars, n_pars, n_groups, iters), dimnames = list(par_names, par_names, groups, NULL)),
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_sig = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    stage = array(stage, iters),
    subj_ll = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    group_a_half = array(NA_real_,dim = c(n_pars, n_groups, iters),dimnames = list(par_names, groups, NULL)),
    a_half = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL))
  )
}

unwind <- function(var_matrix, ...) {
  y <- t(chol(var_matrix))
  diag(y) <- log(diag(y))
  return(y[lower.tri(y, diag = TRUE)])
}


set_epsilon <- function(n_pars, epsilon) {
  if (checkmate::test_null(epsilon)) {
    if (n_pars > 15) {
      epsilon <- 0.1
    } else if (n_pars > 10) {
      epsilon <- 0.3
    } else {
      epsilon <- 0.5
    }
    message(
      sprintf(
        "Epsilon has been set to %.1f based on number of parameters",
        epsilon
      )
    )
  }
  epsilon
}

set_mix <- function(stage, mix) {
  if (checkmate::test_null(mix)) {
    if (stage == "sample") {
      mix <- c(0.1, 0.2, 0.2, 0.7)
    } else {
      mix <- c(0.2, 0.4, 0.4, 0)
    }
    message(
      sprintf(
        "mix has been set to c(%s) based on the stage being run",
        paste(mix, collapse = ", ")
      )
    )
  }
  mix
}

last_sample_group <- function(store, group) {
  list(
    tmu = store$group_mu[,group, store$idx],
    tsig = store$group_sig[, , group,store$idx],
    alpha = store$alpha[, , store$idx],
    tsinv = store$last_group_theta_sig_inv[,,group],
    a_half = store$group_a_half[, group,store$idx]
  )
}

last_sample <- function(store) {
  list(
    tmu = store$theta_mu[, store$idx],
    tsig = store$theta_sig[, , store$idx],
    alpha = store$alpha[, , store$idx],
    tsinv = store$last_theta_sig_inv,
    a_half = store$a_half[, store$idx]
  )
}

numbers_from_proportion <- function(mix_proportion, num_particles = 1000) {
  numbers <- numeric(4)
  draws <- table(sample(c(1:4), size = num_particles, replace = T, prob = mix_proportion))
  numbers[as.numeric(names(draws))] <- draws
  numbers
}

accept_rate <- function(pmwgs, window_size = 200) {
  n_samples <- pmwgs$samples$idx
  if (is.null(n_samples) || n_samples < 3) {
    return(array(0, dim(pmwgs$samples$alpha)[2]))
  }
  if (n_samples <= window_size) {
    start <- 1
    end <- n_samples
  } else {
    start <- n_samples - window_size + 1
    end <- n_samples
  }
  vals <- pmwgs$samples$alpha[1, , start:end]
  apply(
    apply(vals, 1, diff) != 0, # If diff != 0
    2,
    mean
  )
}