#Script for recovering a typical experiment 

rm(list = ls())
library(rtdists)
source("multiGroup/sampling_factors_nGroups_1Cov.R")

log_likelihood=function(x,data, sample=F) {
  x <- exp(x)
  bPars <- grep("b", names(x))
  bs <- x["A"]+x[bPars][data$condition]
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

n.trials = 75       #number trials per subject per conditions
n.subj = 50          #number of subjects
n.cond = 3          #number of conditions


names=c("subject","rt","resp","condition") #names of columns
data = data.frame(matrix(NA, ncol = length(names), nrow = (n.trials*n.subj*n.cond))) #empty data frame
names(data)=names
data$condition = rep(1:n.cond,times = n.trials) #filling in condition
data$subject = rep(1:n.subj, each = n.trials*n.cond) #filling in subjects

parameter.names=c("b1","b2","b3", "A","v1","v2","t0")
n.parameters=length(parameter.names)

ptm1 <- array(dim = n.parameters, dimnames = list(parameter.names)) #an empty array where i will put parameter values
ptm2 <- array(dim = n.parameters, dimnames = list(parameter.names)) #an empty array where i will put parameter values

ptm1[1:n.parameters]=c(0.1,0.3,0.5,0.4,1.2,0.3,-2.4)
ptm2[1:n.parameters]=c(0.1,0.3,0.5,0.4,1.5,0.5,-1.8)
exp(ptm1)
exp(ptm2)
vars = abs(ptm1)/10 #off diagonal correlations are done as absolute/10

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

subj_random_effects_1 <- t(mvtnorm::rmvnorm(n.subj/2,mean=ptm1,sigma=sigmaC))
subj_random_effects_2 <- t(mvtnorm::rmvnorm(n.subj/2,mean=ptm2,sigma=sigmaC))
subj_random_effects <- cbind(subj_random_effects_1, subj_random_effects_2)

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


source("multiGroup/sampling_factors_nGroups_1Cov.R")
data$group <- rep(c(1,2), each = nrow(data)/2)
pars <- rownames(subj_random_effects)

par_idx <- matrix(0, nrow = length(pars), ncol = 2)
par_idx[0:4,] <- "1,2" 
par_idx[5:length(pars),2] <- "1"
par_idx[5:length(pars),1] <- "2"

# Create the Particle Metropolis within Gibbs sampler object ------------------
sampler <- pmwgs(
  data = data,
  pars = pars,
  ll_func = log_likelihood,
  n_factors = 2,
  par_idx = par_idx
)
# start the sampler ---------------------------------------------------------

sampler <- init(sampler, n_cores = 12) # i don't use any start points here

# Sample! -------------------------------------------------------------------
debug(gibbs_step_factor)
burned <- run_stage(sampler, stage = "burn",iter = 500, particles = 100, n_cores =12, pstar = .6)
save(burned, file = paste0("./samples/LBA_", n.subj, "subs_", k, "factors_.RData"))
adapted <- run_stage(burned, stage = "adapt", iter = 5000, particles = 100, n_cores = 12, pstar = .6)
sampled <- run_stage(adapted, stage = "sample", iter = 1000, particles = 100, n_cores = 12, pstar = .6)
save(sampled, file = paste0("./samples/LBA_", n.subj, "subs_", k, "factors_.RData"))

