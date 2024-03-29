#Script for recovering a typical experiment 

rm(list = ls())
source("pmwg/variants/standard.R")

library(rtdists)

log_likelihood=function(x,data, sample=F) {
  x <- exp(x)
  if (sample) { #for sampling
    out=rdiffusion(n=nrow(data),a=x["a"],v=x["v"],t0=x["t0"],z = x["z"],sz = x["sz"], sv = x["sv"], st0 = x["st0"], s = 1)
  } else { #for calculating density
    out=ddiffusion(rt=data$rt,response=data$resp,a=x["a"],v=x["v"],t0=x["t0"],z = x["z"],sz = x["sz"], sv = x["sv"], st0 = x["st0"], s=1)
    bad=(out<1e-10)|(!is.finite(out))
    out[bad]=1e-10
    out=sum(log(out))
  }
  out
}

n.trials = 500       #number trials per subject per conditions
n.subj = 30          #number of subjects
n.cond = 1          #number of conditions


names=c("subject","rt","resp","condition") #names of columns
data = data.frame(matrix(NA, ncol = length(names), nrow = (n.trials*n.subj*n.cond))) #empty data frame
names(data)=names
#data$condition = rep(1:n.cond,times = n.trials) #filling in condition
data$subject = rep(1:n.subj, each = n.trials*n.cond) #filling in subjects

parameter.names=c("a", "v", "t0", "z", "sv", "st0", "sz")
n.parameters=length(parameter.names)

ptm <- array(dim = n.parameters, dimnames = list(parameter.names)) #an empty array where i will put parameter values

ptm[1:n.parameters]=c(0.1, 1, -2, -0.8, 0.5, -2, -1.2)
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

sampler <- pmwgs(
  data = pmwg::forstmann,
  pars = pars,
  ll_func = log_likelihood
)

n_cores = 10

sampler <- init(sampler, n_cores = n_cores) # i don't use any start points here
source("pmwg/variants/standard.R")

burned <- run_stage(sampler, stage = "burn",iter = 500, particles = 100, n_cores = n_cores, pstar = .7)
adapted <- run_stage(burned, stage = "adapt", iter = 1000, particles = 100, n_cores = n_cores, pstar =.7)
debug(variant_funs$get_conditionals)
sampled <- run_stage(adapted, stage = "sample", iter = 500, particles = 100, n_cores = n_cores, pstar = .6)
save(sampled, file = "standard.RData")
