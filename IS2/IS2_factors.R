## Factor IS2 written by Niek, heavily based on code from Reilly Innes and David Gunawan, from Tran et al. 2021 
## set up environment and packages
rm(list=ls())
library(mvtnorm)
library(MCMCpack)
library(rtdists)
library(invgamma)
library(mixtools)
library(condMVNorm)
library(parallel)
library(corpcor) 
library(mvtnorm)
#RJI_change: not sure if this was included before
#library(matrixcalc)

load("~/Documents/pmwg_new/samples/factor_1F_50S_15P_FactRecovery.RData")
data <- sampled$data
sampled$ll_func <- function(pars, data){
  n_pars <- length(pars)
  return(sum(mvtnorm::dmvnorm(as.matrix(data[,2:(1 + n_pars)]), pars, diag(1, nrow = n_pars), log = T)))
}

###### set up variables #####
# number of particles, samples, subjects, random effects etc
n_randeffect=sampled$n_pars
n_subjects = sampled$n_subjects
n_iter = length(sampled$samples$stage[sampled$samples$stage=="sample"])
length_draws = sampled$samples$idx # length of the full transformed random effect vector and/or parameter vector
pars = sampled$par_names

# grab the sampled stage of PMwG
# store the random effects
alpha <- sampled$samples$alpha[,,sampled$samples$stage=="sample"]
# store the mu
theta_mu <- sampled$samples$theta_mu[,sampled$samples$stage=="sample"]
# store the factor loadings. I foresee problems here with the sign switching (if diag unconstrained) and them having normal/t-dist support
lambda <- sampled$samples$theta_lambda[,,sampled$samples$stage=="sample", drop = F]
#store the factor variance matrix and make sure they're on the real line
psi_inv <- sampled$samples$theta_psi_inv[,,sampled$samples$stage=="sample", drop = F]
#store the residual error matrix and make sure they're on the real line
sig_err_inv <- sampled$samples$theta_sig_err_inv[,,sampled$samples$stage=="sample"]
#get the constraints used
constraintMat <- attributes(sampled)$constraintMat

#hyperparameters
n_factors <- ncol(lambda)

unwind_lambda <- function(lambda, constraintMat, n_factors, reverse = F){
  if(reverse){
    out <- matrix(0, n_randeffect, n_factors)
    out[constraintMat] <- lambda
  } else{
    out <- as.numeric(lambda[constraintMat])
  }
  return(out)
}


#unwind lambda
lambda.unwound <- apply(lambda,3,unwind_lambda, constraintMat, n_factors)
sig_err_inv.diag <- log(apply(sig_err_inv, 3, diag))
psi_inv.diag <- matrix(log(apply(psi_inv, 3, diag)), nrow = n_factors)
#Easy to calculate 'faster' but let's be explicit about what we're putting in
n.params<- nrow(alpha) + nrow(theta_mu) + nrow(sig_err_inv.diag) + nrow(psi_inv.diag) + nrow(lambda.unwound)
all_samples=array(dim=c(n_subjects,n.params,n_iter))
mu_tilde=array(dim = c(n_subjects,n.params))
sigma_tilde=array(dim = c(n_subjects,n.params,n.params))

for (j in 1:n_subjects){
  all_samples[j,,] = rbind(alpha[,j,],theta_mu[,],sig_err_inv.diag[,],psi_inv.diag[,],lambda.unwound[,])
  # calculate the mean for re, mu and sigma
  mu_tilde[j,] =apply(all_samples[j,,],1,mean)
  # calculate the covariance matrix for random effects, mu and sigma
  sigma_tilde[j,,] = cov(t(all_samples[j,,]))
}

for(i in 1:n_subjects){ #RJI_change: this bit makes sure that the sigma tilde is pos def
  if(!corpcor::is.positive.definite(sigma_tilde[i,,], tol=1e-8)){
    sigma_tilde[i,,]<-corpcor::make.positive.definite(sigma_tilde[i,,], tol=1e-6)
  }
}

X <- cbind(t(theta_mu),t(sig_err_inv.diag),t(psi_inv.diag), t(lambda.unwound)) 

#perm = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
IS_samples = 250 #number of importance samples
n_particles = 100 #RJI_change: i think 250 should be default

cpus = 8
muX<-apply(X,2,mean)
sigmaX<-var(X)
df=3 #RJI_change: this shouldnt be 1. We should make this an input argument maybe with default at 3 or 5. 1 could be used but its crazy. Also, 30 could be used if wanting a normal dist. 

# generates the 10,000 IS proposals given the mix
prop_theta=mvtnorm::rmvt(IS_samples,sigma = sigmaX, df=df, delta=muX) #RJI_change: uses mix of t dists to generate proposals

group_dist = function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, n_randeffect,
                      constraintMat, n_factors){
  param.theta_mu <- parameters[1:n_randeffect]
  param.sig_err_inv <- exp(parameters[(n_randeffect+1):(n_randeffect + n_randeffect)])
  param.psi_inv <- exp(parameters[(n_randeffect+n_randeffect+1):(n_randeffect + n_randeffect+ n_factors)])
  param.lambda.unwound <- parameters[(n_randeffect+n_randeffect+n_factors+1):length(parameters)] 
  param.lambda <- unwind_lambda(param.lambda.unwound, constraintMat, n_factors, reverse = T)
  param.var <- param.lambda %*% diag(1/param.psi_inv, length(param.psi_inv)) %*% t(param.lambda) + diag(1/param.sig_err_inv)
  if (sample){
    return(mvtnorm::rmvnorm(n_samples, param.theta_mu, param.var))
  }else{
    logw_second<-mvtnorm::dmvnorm(random_effect, param.theta_mu,param.var,log=TRUE)
    return(logw_second)
  }
}

prior_dist = function(parameters, priors = sampled$prior, n_randeffect){ ###mod notes: the sampled$prior needs to be fixed/passed in some other time
  #Extract and when necessary transform back
  param.theta_mu <- parameters[1:n_randeffect]
  param.sig_err_inv <- exp(parameters[(n_randeffect+1):(n_randeffect + n_randeffect)])
  param.psi_inv <- exp(parameters[(n_randeffect+n_randeffect+1):(n_randeffect + n_randeffect+ n_factors)])
  param.lambda.unwound <- parameters[(n_randeffect+n_randeffect+n_factors+1):length(parameters)] 
  #Hyperparameters
  al <- 1 # hyperparameter shape for gamma on psi_inv
  bl <- 1/2 # hyperparameter rate for gamma on psi_inv
  nu <- 2 # hyperparameter shape for gamma on sig_err_inv 
  s2 <- 1/nu # hyperparameter rate for gamma on sig_err_inv
  
  log_prior_mu=sum(dnorm(param.theta_mu, mean = priors$theta_mu_mean, sd = sqrt(priors$theta_mu_var), log =TRUE))
  log_prior_sig_err_inv = sum(dgamma(param.sig_err_inv, shape = 1, rate = .2, log=TRUE))
  log_prior_psi_inv = sum(dgamma(param.psi_inv, shape = 1/2, rate = 1/2, log=TRUE))
  log_prior_lambda=sum(dnorm(param.lambda.unwound, mean = 0, sd = sqrt(priors$theta_lambda_var), log =TRUE))
  
  jac_sig_err_inv <- -sum(log(param.sig_err_inv)) # Jacobian determinant of transformation of log of the sig_err_inv  
  jac_psi_inv <- -sum(log(param.psi_inv)) # Jacobian determinant of transformation of log of the psi_inv
  #Jacobians are actually part of the denominator (dnorm(prop_theta)) since transformations of the data (rather than parameters),
  #Warrant a jacobian added. But we add the jacobians here for ease of calculations. 
  return(log_prior_mu + log_prior_sig_err_inv + log_prior_psi_inv + log_prior_lambda - jac_psi_inv - jac_sig_err_inv)
}

#importance_particles <- array(dim=c(IS_samples,n_particles,n_subjects))

get_logp=function(prop_theta,data,n_subjects,n_particles,n_randeffect,mu_tilde,sigma_tilde,i, group_dist=group_dist,
                  constraintMat, n_factors){
  # make an array for the density
  logp=array(dim=c(n_particles,n_subjects))
  # for each subject, get 1000 IS samples (particles) and find log weight of each
  for (j in 1:n_subjects){
    # generate the particles from the conditional MVnorm AND mix of group level proposals
    wmix <- 0.95
    n1=rbinom(n=1,size=n_particles,prob=wmix)
    if (n1<2) n1=2
    if (n1>(n_particles-2)) n1=n_particles-2 ## These just avoid degenerate arrays.
    n2=n_particles-n1
    # do conditional MVnorm based on the proposal distribution.
    # contrary to standard IS2, all group level pars are 'relevant' and are therefore added to the n.params. We therefore take all prop_theta.
    # for given.ind we only remove the random effects from all parameters. 
    conditional = condMVNorm::condMVN(mean=mu_tilde[j,],sigma=sigma_tilde[j,,],dependent.ind=c(1:n_randeffect),
                                      given.ind=c((n_randeffect+1):n.params),X.given=prop_theta[i,])
    particles1 <- mvtnorm::rmvnorm(n1, conditional$condMean,conditional$condVar)
    # # mix of proposal params and conditional
    particles2 <- group_dist(n_samples=n_particles, parameters = prop_theta[i,],sample=TRUE, n_randeffect=n_randeffect, constraintMat = constraintMat, 
                             n_factors = n_factors)
    particles <- rbind(particles1,particles2)
    for (k in 1:n_particles){
      x <-particles[k,]
      # names for ll function to work
      # mod notes: this is the bit the prior effects
      names(x)<-pars
      # do lba log likelihood with given parameters for each subject, gets density of particle from ll func
      logw_first=sampled$ll_func(x,data = data[as.numeric(factor(data$subject))==j,]) #mod notes: do we pass this in or the whole sampled object????
      # below gets second part of equation 5 numerator ie density under prop_theta
      # particle k and big vector of things
      logw_second <- group_dist(random_effect = particles[k,], parameters = prop_theta[i,], sample= FALSE, n_randeffect = n_randeffect, constraintMat = constraintMat, n_factors = n_factors) #mod notes: group dist
      # below is the denominator - ie mix of density under conditional and density under pro_theta
      logw_third <- log(wmix*mvtnorm::dmvnorm(particles[k,], conditional$condMean, conditional$condVar) + (1-wmix) * exp(logw_second)) #mod notes: fine?
      #logw_third <- log(logw_second) #mod notes: fine?
      # does equation 5
      logw=(logw_first+logw_second)-logw_third
      # assign to correct row/column
      if(is.numeric(logw)){ #RJI_change: protection against weird errors
        logp[k,j]=logw
      }else{
        logp[k,j]=1e-10
      }
      #importance_particles[i,k,j]<-logw
    }
  }
  # we use this part to centre the logw before adding back on at the end. This avoids inf and -inf values
  sub_max = apply(logp,2,max)
  logw = t(t(logp) - sub_max)
  w = exp(logw)
  subj_logp = log(apply(w,2,mean))+sub_max #means
  
  
  # sum the logp and return 
  if(is.nan(sum(subj_logp))){ #RJI_change: this protects against bad particles. Need to check if this should be -Inf or 0 (log scale dependent). i think this is correct though
    return(1e-10)
  }else{
    return(sum(subj_logp))
  }
}

compute_lw=function(prop_theta,data,n_subjects,n_particles,n_randeffect,mu_tilde,sigma_tilde,i, prior_dist=prior_dist, sampled=sampled, constraintMat, n_factors){
  
  logp.out <- get_logp(prop_theta, data, n_subjects, n_particles, n_randeffect, mu_tilde, sigma_tilde, i, group_dist = group_dist, constraintMat, n_factors)
  #do equation 10
  logw_num <- logp.out+prior_dist(parameters = prop_theta[i,], priors = sampled$prior, n_randeffect)
  logw_den <- mvtnorm::dmvt(prop_theta[i,], delta=muX, sigma=sigmaX,df=df, log = TRUE) #RJI_change: no longer the mix, just the density of multi t
  logw <- logw_num-logw_den # this is the equation 10
  return(c(logw))
  ##NOTE: we should leave a note if variance is shit - variance is given by the logp function (currently commented out)
}

##### make it work

#makes an array to store the IS samples
tmp<-array(dim=c(IS_samples))

#do the sampling
if (cpus>1){
  tmp <- mclapply(X=1:IS_samples,mc.cores = cpus, FUN = compute_lw, prop_theta = prop_theta,data = data,n_subjects= n_subjects,n_particles = n_particles,
                  n_randeffect = n_randeffect,mu_tilde=mu_tilde,sigma_tilde = sigma_tilde, prior_dist=prior_dist, sampled=sampled, constraintMat, n_factors)
} else{
  for (i in 1:IS_samples){
    cat(i)
    # save.image("tmpbug10.RData")
    tmp[i]<-compute_lw(prop_theta, data, n_subjects, n_particles, n_randeffect, mu_tilde, sigma_tilde, i, prior_dist=prior_dist, sampled=sampled, constraintMat, n_factors)
  }
}

finished <- tmp
tmp<-unlist(tmp)
max.lw <- max(tmp)
mean.centred.lw <- mean(exp(tmp-max.lw)) #takes off the max and gets mean (avoids infs)
lw <- log(mean.centred.lw)+max.lw #puts max back on to get the lw
lw



