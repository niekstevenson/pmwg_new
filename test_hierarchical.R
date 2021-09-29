rm(list = ls())
library(rtdists)
source("pmwg/sampling.R")

log_likelihood=function(x,data,sample=TRUE) {
  x=exp(x)
  bs=x["A"]+x[c("b1","b2","b3")][data$condition]
  if (sample) { #for sampling
    out=rLBA(n=nrow(data),A=x["A"],b=bs,t0=x["t0"],mean_v=x[c("v1","v2")],sd_v=c(1,1),distribution="norm",silent=TRUE)
  } else { #for calculating density
    out=dLBA(rt=data$rt,response=data$resp,A=x["A"],b=bs,t0=x["t0"],mean_v=x[c("v1","v2")],sd_v=c(1,1),distribution="norm",silent=TRUE)
    bad=(out<1e-10)|(!is.finite(out))
    out[bad]=1e-10
    out=sum(log(out))
  }
  out
}

n.trials = 75      #number trials per subject per conditions
n.subj = 5        #number of subjects
n.cond = 3          #number of conditions


names=c("subject","rt","resp","condition") #names of columns
data = data.frame(matrix(NA, ncol = length(names), nrow = (n.trials*n.subj*n.cond))) #empty data frame
names(data)=names
data$condition = rep(1:n.cond,times = n.trials) #filling in condition
data$subject = rep(1:n.subj, each = n.trials*n.cond) #filling in subjects

parameter.names=c("b1","b2","b3", "A","v1","v2","t0")
n.parameters=length(parameter.names)

ptm <- array(dim = n.parameters, dimnames = list(parameter.names)) #an empty array where i will put parameter values

ptm[1:n.parameters]=c(0.1,0.3,0.5,0.4,1.2,0.3,-2.4)
exp(ptm)
vars = abs(ptm)/10 #off diagonal correlations are done as absolute/10

sigmaC <- matrix(c(.8, .5, .4, .15, .15, .3, -.15,
                   .5, .8, .4, .2, .3, .3, .3,
                   .4, .4, .8, .1, .1, .2, .2, 
                   .15, .2, .1, .8, .2, .2, .1, 
                   .15, .3, .1, .2, .8, .5, .2, 
                   .3, .3, .2, .2, .5, .8, .1,
                   -.15, .3, .2, .1, .2, .1, .8),
                 nrow=7,ncol=7)

###std dev correlation on diagonal - you might think this should be corr = 1, but it's actually the standard deviation 
diag(sigmaC)=sqrt(vars)
sigmaC <- sdcor2cov(sigmaC)

subj_random_effects <- t(mvtnorm::rmvnorm(n.subj,mean=ptm,sigma=sigmaC))


for (i in 1:n.subj){
  tmp<- log_likelihood(subj_random_effects[,i],sample=TRUE,data=data[data$subject==i,])
  data$rt[data$subject==i]=tmp$rt
  data$resp[data$subject==i]=tmp$response
}

pars <- rownames(subj_random_effects)

priors <- list(
  theta_mu_mean = rep(0, length(pars)),
  theta_mu_var = diag(rep(1, length(pars)))
)


# Create the Particle Metropolis within Gibbs sampler object ------------------
lba_loglike <- function(x, data, sample = FALSE) {
  x <- exp(x)
  if (any(data$rt < x["t0"])) {
    return(-1e10)
  }
  if (sample){
    data$rt=NA
    data$resp = NA
  }
  
  bs <- x["A"] + x[c("b1", "b2", "b3")][data$condition]
  
  if (sample) {
    out <- rtdists::rLBA(n = nrow(data),
                         A = x["A"],
                         b = bs,
                         t0 = x["t0"],
                         mean_v = x[c("v1", "v2")],
                         sd_v = c(1, 1),
                         distribution = "norm",
                         silent = TRUE)
    data$rt <- out$rt
    data$resp <- out$resp
    
  } else {
    out <- rtdists::dLBA(rt = data$rt,
                         response = data$resp,
                         A = x["A"],
                         b = bs,
                         t0 = x["t0"],
                         mean_v = list(x["v1"],x[ "v2"]),
                         sd_v = c(1, 1),
                         distribution = "norm",
                         silent = TRUE)
    bad <- (out < 1e-10) | (!is.finite(out))
    out[bad] <- 1e-10
    out <- sum(log(out))
  }
  if (sample){return(data)}
  if (!sample){return(out)}
}

sampler <- pmwgs(
  data = data,
  pars = pars,
  prior = priors,
  ll_func = lba_loglike
)
sampler <- init(sampler, n_cores = 1) # i don't use any start points here
# Sample! -------------------------------------------------------------------
burned <- run_stage(sampler, stage = "burn",iter = 1000, particles = 100, n_cores = 16, pstar = .7)
adapted <- run_stage(burned, stage = "adapt", iter = 5000, particles = 100, n_cores = 16, pstar =.6)
sampled <- run_stage(adapted, stage = "sample", iter = 2500, particles = 100, n_cores = 16, pstar = .6)