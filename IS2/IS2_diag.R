## diag IS2 written by Niek Stevenson heavily inspired by Reilly Innes and David Gunawan, from Tran et al. 2021 
## set up environment and packages 
rm(list=ls())
library(mvtnorm)
library(MCMCpack)
library(rtdists)
library(invgamma)
library(mixtools)
library(condMVNorm)
library(parallel)
library(corpcor) #RJI_change: not sure if this was included before
#library(matrixcalc)

load("~/Documents/UVA/2021/pmwg_new/samples/LBA_25subs_0factors_.RData")
data <- sampled$data

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
theta <- sampled$samples$theta_mu[,sampled$samples$stage=="sample"]
# store the cholesky transformed sigma
sig <- sampled$samples$theta_var[,,sampled$samples$stage=="sample"]
# the a-half is used in  calculating the Huang and Wand (2013) prior. 
# The a is a random sample from inv gamma which weights the inv wishart. The mix of inverse wisharts is the prior on the correlation matrix
a_half <- log(sampled$samples$a_half[,sampled$samples$stage=="sample"])

unwind <- function(x,reverse=FALSE, diag = TRUE) {
  if (reverse) {
    if(diag){
      out <- diag(exp(x), nrow = lenght(x))
    } else{
      out <- exp(x)
    }
  } else {
    out <- log(diag(x))
  }
  return(out)
}

#unwound sigma
pts2.unwound = apply(sig,3,unwind)
n.params<- nrow(theta) + nrow(pts2.unwound)+ nrow(a_half)
all_samples=array(dim=c(n_subjects,n.params,n_iter))
mu_tilde=array(dim = c(n_subjects,n.params))
sigma_tilde=array(dim = c(n_subjects,n.params,n.params))

for (j in 1:n_subjects){
  all_samples[j,,] = rbind(alpha[,j,],theta[,],pts2.unwound[,])
  # calculate the mean for re, mu and sigma
  mu_tilde[j,] =apply(all_samples[j,,],1,mean)
  # calculate the covariance matrix for random effects, mu and sigma
  sigma_tilde[j,,] = cov(t(all_samples[j,,]))
}

X <- cbind(t(theta),t(pts2.unwound),t(a_half)) 

#perm = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
IS_samples = 250 #number of importance samples
n_particles = 100 #RJI_change: i think 250 should be default

cpus = 12
muX<-apply(X,2,mean)
sigmaX<-var(X)
df=3 #RJI_change: this shouldnt be 1. We should make this an input argument maybe with default at 3 or 5. 1 could be used but its crazy. Also, 30 could be used if wanting a normal dist. 

# generates the 10,000 IS proposals given the mix
prop_theta=mvtnorm::rmvt(IS_samples,sigma = sigmaX, df=df, delta=muX) #RJI_change: uses mix of t dists to generate proposals

group_dist = function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, n_randeffect){
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)] 
  param.theta.sig2 <- unwind(param.theta.sig.unwound, reverse = TRUE) 
  #dmvnorm will be slow, because since it's a diagonal matrix, multiple dnorm calls would be quicker
  if (sample){
    return(mvtnorm::rmvnorm(n_samples, param.theta.mu,param.theta.sig2))
  }else{
    logw_second<-mvtnorm::dmvnorm(random_effect, param.theta.mu,param.theta.sig2,log=TRUE)
    return(logw_second)
  }
}

prior_dist = function(parameters, prior_parameters = sampled$prior, n_randeffect){ ###mod notes: the sampled$prior needs to be fixed/passed in some other time
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)] 
  param.theta.sig2 <- unwind(param.theta.sig.unwound, reverse = TRUE, diag = FALSE)
  param.a <- exp(parameters[((length(parameters)-n_randeffect)+1):(length(parameters))])
  v_half <- 2
  a_half <- 1
  log_prior_mu=mvtnorm::dmvnorm(param.theta.mu, mean = prior_parameters$theta_mu_mean, sigma = prior_parameters$theta_mu_var, log =TRUE)
  log_prior_sigma = sum(invgamma::dinvgamma(param.theta.sig2, shape = v_alpha/2, rate = v_alpha/param.a, log = T))
  log_prior_a = sum(invgamma::dinvgamma(param.a,scale = 1/2,shape=1/(a_half^2),log=TRUE))
  jac_a <- -sum(log(param.a))
  jac_sig <- -sum(log(1/param.theta.sig2))
  return(log_prior_mu + log_prior_sigma + log_prior_a - jac_a - jac_sig)
}

get_logp=function(prop_theta,data,n_subjects,n_particles,n_randeffect,mu_tilde,sigma_tilde,i, group_dist=group_dist){
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
    # do conditional MVnorm based on the proposal distribution
    conditional = condMVNorm::condMVN(mean=mu_tilde[j,],sigma=sigma_tilde[j,,],dependent.ind=c(1:n_randeffect),
                                      given.ind=c((n_randeffect+1):n.params),X.given=prop_theta[i,c(1:(n.params-n_randeffect))])
    particles1 <- mvtnorm::rmvnorm(n1, conditional$condMean,conditional$condVar)
    # # mix of proposal params and conditional
    particles2 <- group_dist(n_samples=n_particles, parameters = prop_theta[i,],sample=TRUE, n_randeffect=n_randeffect)
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
      logw_second <- group_dist(random_effect = particles[k,], parameters = prop_theta[i,], sample= FALSE, n_randeffect = n_randeffect) #mod notes: group dist
      # below is the denominator - ie mix of density under conditional and density under pro_theta
      logw_third <- log(wmix*dmvnorm(particles[k,], conditional$condMean, conditional$condVar) + (1-wmix) * exp(logw_second)) #mod notes: fine?
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

compute_lw=function(prop_theta,data,n_subjects,n_particles,n_randeffect,mu_tilde,sigma_tilde,i, prior_dist=prior_dist, sampled=sampled){
  
  logp.out <- get_logp(prop_theta, data, n_subjects, n_particles, n_randeffect, mu_tilde, sigma_tilde, i, group_dist = group_dist)
  #do equation 10
  logw_num <- logp.out[1]+prior_dist(parameters = prop_theta[i,], prior_parameters = sampled$prior, n_randeffect)
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
                  n_randeffect = n_randeffect,mu_tilde=mu_tilde,sigma_tilde = sigma_tilde, prior_dist=prior_dist, sampled=sampled)
} else{
  for (i in 1:IS_samples){
    cat(i)
    # save.image("tmpbug10.RData")
    tmp[i]<-compute_lw(prop_theta, data, n_subjects, n_particles, n_randeffect, mu_tilde, sigma_tilde, i, prior_dist=prior_dist, sampled=sampled)
  }
}


finished <- tmp
tmp<-unlist(tmp)
max.lw <- max(tmp)
mean.centred.lw <- mean(exp(tmp-max.lw)) #takes off the max and gets mean (avoids infs)
lw <- log(mean.centred.lw)+max.lw #puts max back on to get the lw
lw



