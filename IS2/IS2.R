## standard IS2 written by Reilly Innes and David Gunawan, from Tran et al. 2021 
## set up environment and packages 
library(mvtnorm)
library(MCMCpack)
library(rtdists)
library(invgamma)
library(mixtools)
library(condMVNorm)
library(parallel)
library(corpcor) #RJI_change: not sure if this was included before
#library(matrixcalc)


IS2 <- function(samples, filter, IS_samples = 1000, n_particles = 250, n_cores = 1, df = 3){
  ###### set up variables #####
  info <- add_info_base(samples)
  all_pars <- variant_funs$get_all_pars(samples, filter, info)
  muX<-apply(all_pars$X,2,mean)
  varX<-cov(all_pars$X)

  prop_theta=mvtnorm::rmvt(IS_samples,sigma = varX, df=df, delta=muX)

  #do the sampling
  logw_num <- mclapply(X=1:IS_samples, 
                  FUN = compute_lw_num, 
                  prop_theta = prop_theta,
                  n_particles = n_particles,
                  mu_tilde=all_pars$mu_tilde,
                  var_tilde = all_pars$var_tilde, 
                  info = all_pars$info,
                  mc.cores = n_cores)
  logw_den <- mvtnorm::dmvt(prop_theta, delta=muX, sigma=varX,df=df, log = TRUE)
  
  finished <- unlist(logw_num) - logw_den
  max.lw <- max(finished)
  mean.centred.lw <- mean(exp(finished-max.lw)) #takes off the max and gets mean (avoids infs)
  lw <- log(mean.centred.lw)+max.lw #puts max back on to get the lw
  lw
}

get_logp=function(prop_theta,n_particles,mu_tilde,var_tilde, info){
  # Unload for legibility
  n_subjects <- info$n_subjects
  n_randeffect <- info$n_randeffect
  n_params <- info$n_params
  ll_func <- info$ll_func
  data <- info$data
  # make an array for the density
  logp <- array(dim=c(n_particles,n_subjects))

  # for each subject, get 1000 IS samples (particles) and find log weight of each
  for (j in 1:n_subjects){
    # generate the particles from the conditional MVnorm AND mix of group level proposals
    wmix <- 0.95
    n1=rbinom(n=1,size=n_particles,prob=wmix)
    if (n1<2) n1=2
    if (n1>(n_particles-2)) n1=n_particles-2 ## These just avoid degenerate arrays.
    n2=n_particles-n1
    # do conditional MVnorm based on the proposal distribution
    conditional = condMVNorm::condMVN(mean=mu_tilde[j,],sigma=var_tilde[j,,],
                                      dependent.ind=c(1:n_randeffect),
                                      given.ind=c((n_randeffect+1):n_params),
                                      X.given=prop_theta[c(1:(n_params-n_randeffect))],
                                      check.sigma = F)
    # Subject specific efficient distribution
    particles1 <- mvtnorm::rmvnorm(n1, conditional$condMean,conditional$condVar)
    # Group level
    particles2 <- variant_funs$group_dist(n_samples=n2, parameters = prop_theta,
                             sample=TRUE, info = info)
    particles <- rbind(particles1,particles2)
    # names for ll function to work
    colnames(particles) <- info$par_names
    # do lba log likelihood with given parameters for each subject, gets density of particle from ll func
    lw_first <- apply(particles, 1, ll_func, data = data[data$subject==unique(data$subject)[j],])
    # below gets second part of equation 5 numerator ie density under prop_theta
    lw_second <- apply(particles, 1, variant_funs$group_dist, prop_theta, FALSE, NULL, info)
    # below is the denominator - ie mix of density under conditional and density under pro_theta
    lw_third <- log(wmix*dmvnorm(particles, conditional$condMean, conditional$condVar) + (1-wmix) * exp(lw_second)) 
    # does equation 5
    logw=lw_first+lw_second-lw_third
    logw[!is.numeric(logw)] <- 1e-10
    logp[,j]=logw
  }
  # we use this part to centre the logw before adding back on at the end. This avoids inf and -inf values
  # Niek has some doubts about whether this is correct
  sub_max = apply(logp,2,max)
  logw = sweep(logp, 2, sub_max)
  w = exp(logw)
  subj_logp = log(apply(w,2,mean))+sub_max #means
  
  # sum the logp and return 
  if(is.nan(sum(subj_logp))){ 
    return(1e-10)
  }else{
    return(sum(subj_logp))
  }
}

compute_lw_num=function(i, prop_theta,n_particles,mu_tilde,var_tilde,info){
  logp.out <- get_logp(prop_theta[i,], n_particles, mu_tilde, var_tilde, info)
  logw_num <- logp.out+variant_funs$prior_dist(parameters = prop_theta[i,], info)
  return(logw_num)
}

add_info_base <- function(samples){
  info <- list(
    n_randeffect = samples$n_pars,
    n_subjects = samples$n_subjects,
    data = samples$data,
    ll_func = samples$ll_func,
    prior = samples$prior,
    hyper = attributes(samples)
  )
  return(info)
}

