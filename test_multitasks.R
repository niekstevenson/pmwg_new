#Script for testing the recovery of joint models consisting of multiple tasks (can be done in standard or factors or whatever)

rm(list = ls())
source("pmwg/variants/standard.R")
library(rtdists)

joint_ll <- function(x, data){
  parPreFixs <- gsub("[|].*", "", names(x))
  totalSum <- 0
  for(mod in unique(parPreFixs)){
    currentPars <- x[which(parPreFixs == mod)]
    names(currentPars) <- gsub(".*[|]", "", names(currentPars))
    modelData <- data[[mod]][[1]]
    totalSum <- totalSum + log_likelihood(currentPars, modelData, sample = F)
  }
  return(totalSum)
}

log_likelihood=function(x,data, sample=TRUE) {
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


n.trials <- 100      #number trials per subject per conditions
n.subj <- 25 #number of subjects
n.cond <- 3
n.exp <- 2

allparameters <- numeric()
alldata <- list()


for(i in 1:n.exp){
  names=c("subject","rt","resp","condition") #names of columns
  data <- data.frame(matrix(NA, ncol = length(names), nrow = (n.trials*n.subj*n.cond))) #empty data frame
  names(data) <- names
  data$condition <- rep(1:n.cond,times = n.trials) #filling in condition
  data$subject <- rep(1:n.subj, each = n.trials*n.cond) #filling in subjects
  alldata[[i]] <- data
  parameter.names <- c(c("A","v1","v2","t0"), c(paste0("b", 1:n.cond)))
  parameters <- c((0.1 + .02*i)*(-1)^i + c(0.4, 1.2, .7, -2.4), seq(0.1, 0.5, length.out = n.cond))
  names(parameters) <- paste0("Mod", i, "|", parameter.names)
  allparameters <- c(allparameters, parameters)
}

fillCorMatrix <- function(matrix, name1, name2, value){
  matrix[grep(name1, rownames(matrix)), grep(name2, colnames(matrix))] <- value
  matrix[grep(name2, rownames(matrix)), grep(name1, colnames(matrix))] <- value
  return(matrix)
}

adjustCorMatrix <- function(matrix, adj){
  #Was probably a cleverer way to do this
  rowNames <- rownames(matrix)
  for(i in 1:nrow(matrix)){
    rowPreFix <- substr(rowNames[i], 0, 4)
    matrix[i,!grepl(rowPreFix, colnames(matrix))] <- .3*matrix[i,!grepl(rowPreFix, colnames(matrix))] 
  }
  return(matrix)
}

exp(allparameters)
n.parameters=length(allparameters)
corMat <- matrix(0, nrow = n.parameters, ncol = n.parameters)
rownames(corMat) <- colnames(corMat) <- names(allparameters)
corMat[grep("A", rownames(corMat)),] <- .3 
corMat[,grep("A", colnames(corMat))] <- .3 

corMat <- fillCorMatrix(corMat, "v", "b", .25)
corMat <- fillCorMatrix(corMat, "v", "v", .4)
corMat <- fillCorMatrix(corMat, "b", "b", .4)
corMat <- fillCorMatrix(corMat, "t0", "v", -.1)
corMat <- fillCorMatrix(corMat, "t0", "b", -.25)
corMat <- adjustCorMatrix(corMat, .33)



vars = abs(allparameters)/10 #off diagonal correlations are done as absolute/10

###std dev correlation on diagonal - you might think this should be corr = 1, but it's actually the standard deviation 
diag(corMat)=sqrt(vars)
corMat <- apply(burned$samples$theta_var, 1:2, mean) #sdcor2cov(corMat)

subj_random_effects <- t(mvtnorm::rmvnorm(n.subj,mean=allparameters,sigma=corMat))
exp(subj_random_effects)
parPreFixs <- gsub("[|].*", "", rownames(subj_random_effects))
df <- data.frame(subject = 1:n.subj)
for(i in 1:n.exp){
  modIdx <- which(parPreFixs == paste0("Mod", i))
  task_random_effects <- subj_random_effects[modIdx,]
  rownames(task_random_effects) <- gsub(".*[|]", "", rownames(task_random_effects))
  data <- alldata[[i]]
  for (j in 1:n.subj){
    tmp<- log_likelihood(task_random_effects[,j],sample=TRUE,data=data[data$subject==j,])
    data$rt[data$subject==j]=tmp$rt
    data$resp[data$subject==j]=tmp$response
  }
  df[paste0("Mod", i)] <- I(list(split(data, f = data$subject)))
}


pars <- rownames(subj_random_effects)
# Create the Particle Metropolis within Gibbs sampler object ------------------
sampler <- pmwgs(
  data = df,
  pars = pars,
  ll_func = joint_ll,
)
# start the sampler ---------------------------------------------------------
sampler <- init(sampler, n_cores = 12) # i don't use any start points here

# Sample! -------------------------------------------------------------------
burned <- run_stage(sampler, stage = "burn",iter = 2000, particles = 100, n_cores =12, pstar = .6)
adapted <- run_stage(burned, stage = "adapt", iter = 5000, particles = 100, n_cores = 12, pstar = .6)
sampled <- run_stage(adapted, stage = "sample", iter = 1000, particles = 100, n_cores = 12, pstar = .6)

