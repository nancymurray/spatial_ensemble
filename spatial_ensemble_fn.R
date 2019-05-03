########################################################################################
#
#         Spatial Ensemble Averaging with Statistical Models
#             
#
#                   Nancy L. Murray
#                 nancy.murray@emory.edu
#
#                   April 29, 2019
#
# Reference: 
# "A Bayesian Ensemble Approach to Combine PM2.5 Estimates from Statistical Models using Satellite
#Imagery and Numerical Model Simulation"
#Nancy L Murray, Heather A Holmes, Yang Liu, Howard H Chang
#
#
#
# MODEL: 
#### y = p*m1 + p*m2
#### W = logit(p)
#### W ~ N (0, S), S[i,j] = tau2 * exp (-1/rho * d_{ij})
#### i,j = 1, 2, ..., S (number of spatial locations)

# PRIORS:
# 
#### tau2 ~ Inv-gamma (a, b); 
#### rho ~ gamma (c, d)
#

########################################################

##################################################################################################
# Function to obtain spatial weights (spatial_ensemble)
##################################################################################################

###
###############
### INPUTS  ###
###############
### DATA

# train.cmaq.dat - training data that must have: 
#      (1) mean outputs from Bayesian hierarchical models (BHM) for CMAQ "Pred"
#      (2) standard error outputs from Bayesian hierarchical models (BHM) for CMAQ "Pred_SD"
#      (3) observed PM values as "Obs"
#      (4) monitor IDs as "MonID"
#train.aod.dat - training data that must have: 
  #      (1) mean outputs from Bayesian hierarchical models (BHM) for AOD "Pred"
  #      (2) standard error outputs from Bayesian hierarchical models (BHM) for AOD "Pred_SD"
  #      (3) observed PM values as "Obs"
  #      (4) monitor IDs as "MonID"


### SPACE-TIME INDEX
### Dist.mat = S x S matrix of pairwise Euclidean distance where S is the total number of locations


# DEFAULT PRIORS:
# These are chosen to be uninformative. 
#### tau2 ~ Inv-gamma (0.001, 0.001); #can change with a and b
#### rho ~ gamma (5, 0.05) #can change with c and d


# Other Model specifications
#tau2= 1; starting value for tau2
# tunevar=0.45; tuning 
# rho = 10; starting value for rho
# rho.tune=0.4; tuning


### MCMC ITERATIONS
### n.iter = number of "total" MCMC iterations
### n.burn = number of initial burn "discarded" interations


###############
### OUTPUT  ###
###############
###
#logit_weight: logit weights saved for each location s
#ensemble_weight: ensemble weights saved for each location s
#tau2: posterior samples of tau
#rho: posterior samples of rho


n.iter = 10000 #number of MCMC iterations to run
n.burn = 2500 #number of burn in 
postburn <- n.iter - n.burn 
rho.tune=0.4
tunevar=0.45


#Initialize
tau2= 1
rho = 10

#Setting the tau and rho limits
a = 0.001
b = 0.001
c = 5
d = 0.05



########################################################
########################################################
########################################################
# Set up an empty results dataframe
JointFullResults <- NULL

# Set a seed for random number generation
set.seed(0428)


spatial_ensemble <- function(train.cmaq.dat=train.cmaq.dat, train.aod.dat=train.aod.dat, Dist.mat=Dist.mat,
n.iter = 10000, 
n.burn = 2500, 
postburn = n.iter - n.burn, 
rho.tune=0.4, tunevar=0.45, tau2= 1,
rho = 10,
a = 0.001,
b = 0.001,
c = 5,
d = 0.05){

N.loc = nrow(Dist.mat)


#Save posterior samples
tau2.save = rho.save  = rep (NA,(n.iter+1))
tau2.save[1] <- tau2
rho.save[1] <- rho

  
n = nrow (Dist.mat)

#set up saving weight posteriors
p <- rep(0.5, N.loc)
p.save <- matrix(NA, ncol=(n.iter+1), nrow=N.loc)
p.save[, 1] <- p
rho.acc = 0
p.acc = rep(0, N.loc)
#TRANSFORM THE p TO BE APPROXIMATELY NORMAL
logitp <- log(p/(1-p))
W <- logitp
W.save <- matrix(NA, ncol=(n.iter+1), nrow=N.loc)
W.save[,1] <- W


R = exp(-1/rho*Dist.mat)

monlocvec <- sort(unique(train.cmaq.dat$MonID))


S <- nrow(Dist.mat)
##################################################################################################
##Begin MCMC for CV
#Update the weights at each iteration and use results from iterations to estimate PM.
##################################################################################################

for (i in 1:n.iter){
  #i = 1
  if (i %% 1000 == 0){ print (i)}
  
  
  big.sig <- tau2 * exp (-1/rho * Dist.mat)
  
  

  for (s in 1:S) { #begin training locations loop
    
    
    ResultsCMAQ <- train.cmaq.dat[which(train.cmaq.dat$MonID==monlocvec[s]),]
    ResultsAOD <- train.aod.dat[which(train.aod.dat$MonID==monlocvec[s]),]
    
    if (i %% 1000 == 0){ print(paste("Location", s, "of", S)) }
    ptm <- proc.time()
    
    
    #Specific to each location
    days <- length(ResultsCMAQ$Pred)
    zt <- matrix(0,nrow=(n.iter+1), ncol=days)
    
    
    
    #initialize z and weight p
    zt[1,] <- 1
    
    #means and standard errors
    m1 <- ResultsCMAQ$Pred
    m2 <- ResultsAOD$Pred
    se1 <- ResultsCMAQ$Pred_SD
    se2 <- ResultsAOD$Pred_SD 
    y <- ResultsCMAQ$Obs
    
    d1 <- c()
    d2 <- c()
    pt <- c()
    
    
    for (k in 1:days) { #This is the loop for T days at each s
      d1[k] <- dnorm(y[k], m1[k], se1[k])
      d2[k] <- dnorm(y[k], m2[k], se2[k])
      pt[k] <- d1[k]*p[s]/(d1[k]*p[s] + d2[k]*(1-p[s]))
      
      zt[i+1,k] <- rbinom(1,1,pt[k])
      
    } # end of k in 1 through T days loop 
    
    #Update W (W is just the logit of the weight here)
    
    W_s <- W[s]
    
    
    #Conditionals
    sig11 <- big.sig[s,s]
    sig12 <- big.sig[s,-s]
    sig21 <- big.sig[-s,s]
    sig22inv <- solve(big.sig[-s,-s]) #ginv(big.sig[-s,-s])
    wvec <- W[-s] #W is already in logit form
    condmean <- sig12%*%sig22inv%*%wvec
    condvar <- sig11 - sig12%*%sig22inv%*%sig21
    
    
    #propose new weight
    wt.prop = rnorm(1, W_s, sqrt(tunevar))
    mu.tilde = condmean
    sig.tilde = sqrt(condvar)
    
    
    #transforming the first part (binomial) to be in logit form
    lik.prop = sum(zt[i+1,]*wt.prop) - days*log(1+exp(wt.prop)) + dnorm(wt.prop, mu.tilde, sig.tilde, log=T) 
    lik.curr = sum(zt[i+1,]*W_s) - days*log(1+exp(W_s)) + dnorm(W_s,mu.tilde, sig.tilde, log=T)
    
    acc.prob = min (0, lik.prop-lik.curr, na.rm=T)
    
    if ( log(runif(1)) < acc.prob){
      p[s] = expit(wt.prop)
      
      p.acc[s] = p.acc[s] + 1/n.iter 
      
      W[s] <- wt.prop
    }
    
    
    W.save[s,i+1] = W[s]
    p.save[s,i+1] = p[s]

    
  } # end training locations loop; end loop of s locations
  
  
  
  #Update tau2
  SS = t(W)%*%solve(R)%*%W 
  tau2 = rinvgamma (1, a + n/2, b + 1/2*SS)
  
  #Update rho
  rho.prop = rlnorm(1, log(rho), rho.tune)
  R.prop = exp (-(1/rho.prop*Dist.mat) ) 
  lik.prop = -0.5*log(det(R.prop)) - 0.5/tau2*t(W)%*%solve(R.prop)%*%W + (c-1)*log(rho.prop) - d*rho.prop + log(rho.prop)
  lik.curr = -0.5*log(det(R)) - 0.5/tau2*t(W)%*%solve(R)%*%W + (c-1)*log(rho)- d*rho + log(rho)
  acc.prob = min (0, lik.prop-lik.curr, na.rm=T)
  
  if ( log(runif(1)) < acc.prob){
    rho = rho.prop
    rho.acc = rho.acc + 1/n.iter 
    R = R.prop
  } # end evaluation of likelihood loop for rho
  #end update of rho
  
  #Save samples
  tau2.save[i+1] = tau2
  rho.save[i+1] = rho
  
  
  
} #end n.iter loop
print ("MCMC complete")
#####################################################################
print ("Saving post-burn iterations")
tau2.save = tau2.save[-c(1:(n.burn+1))]
rho.save = rho.save[-c(1:(n.burn+1))]
p.save = p.save[,-c(1:(n.burn+1))]
W.save = W.save[,-c(1:(n.burn+1))]


list(logit_weight = W.save, ensemble_weight = p.save, tau2 =tau2.save, rho=rho.save)

} #end spatial ensemble function 


########################################################################################
## Function to perform predictions (pred_spatial_ensemble)
########################################################################################

################
###  INPUTS  ###
################
#Data
# test.cmaq.dat - test data that must have: 
#      (1) mean outputs from Bayesian hierarchical models (BHM) for CMAQ "Pred"
#      (2) standard error outputs from Bayesian hierarchical models (BHM) for CMAQ "Pred_SD"
#      (3) monitor IDs as "MonID"
#test.aod.dat - test data that must have: 
#      (1) mean outputs from Bayesian hierarchical models (BHM) for AOD "Pred"
#      (2) standard error outputs from Bayesian hierarchical models (BHM) for AOD "Pred_SD"
#      (3) monitor IDs as "MonID"
#
# obj = fitted MCMC object with outputs from spatial_ensemble
#
# Dist.mat.new =  S x S matrix of pairwise Euclidean distance where S is the total number of locations
#
#Dist.mat.new = matrix(NA, ncol=N.loc, nrow=N.loc.new)
#
# n.iter = number Monte Carlo realizations
#
#
################
### OUTPUTS  ###
################
#
#
#### spatial_weight: weight for spatial location (based on conditional mean) 
### est_mean: ensemble estimate of the PM2.5 levels for each prediction record
### est_se: standard error of the ensemble estimate


pred_spatial_ensemble <- function(obj, Dist.mat=Dist.mat, Dist.mat.new = Dist.mat.new, test.cmaq.dat=test.cmaq.dat, test.aod.dat=test.aod.dat) {
JointFullResults <- NULL
N.loc.new = nrow(Dist.mat.new) #length of new locations to predict 
N.loc = nrow(Dist.mat) #length of old locations
monlocvectest <- sort(unique(test.cmaq.dat)) #location IDs for predicted locations

sig12 <- matrix(NA, ncol=N.loc, nrow=N.loc.new)
sig22inv <- matrix(NA, ncol=N.loc, nrow=N.loc)
wvec <- c()
condmean <- matrix(NA, ncol=postburn, nrow=N.loc.new)

#assign saved values to object in pred function
tau2.save <- obj$tau2
rho.save <- obj$rho
W.save <- obj$logit_weight

#conditional posterior
for (j in 1:postburn) {
  sig12 <- tau2.save[j] * exp (-1/rho.save[j] * Dist.mat.new)
  sig22inv <- solve(tau2.save[j] * exp (-1/rho.save[j] * Dist.mat))
  wvec <- W.save[,j]
  condmean[,j] <- sig12%*%sig22inv%*%wvec
} #close j in 1:postburn loop

savedparsjoint.m <- NULL
for (s in 1:N.loc.new) { #begin prediction location loop
  savedparsjoint.m <- cbind(condmean[s,], tau2.save, rho.save)
  Whead <- paste("W saved", s)
  colnames(savedparsjoint.m) <- c(Whead, "Tau2","Rho")
  
  
  #need the prediction cmaq and aod for full prediction (after weights obtained)
  ### Remember we need to convert logit weights to weights
  weightmat <- exp(condmean[s,])/(1+exp(condmean[s,]))
  
  postpredmed <- median(weightmat)
  postpredmean <- mean(weightmat)
  meancmaqp <-  mean(expit(condmean[s,])) #proper conversion from logit normal; mean of the entire backtransformed sample
  
  
  ResultsCMAQ <- test.cmaq.dat[which(test.cmaq.dat$MonID==monlocvectest[s]),]
  ResultsAOD <- test.aod.dat[which(test.aod.dat$MonID==monlocvectest[s]),]
  
  
  combinedmu <- meancmaqp*ResultsCMAQ$Pred + (1-meancmaqp)*ResultsAOD$Pred
  combinedse <- ((meancmaqp^2)*(ResultsCMAQ$Pred_SD^2) + ((1-meancmaqp)^2)*(ResultsAOD$Pred_SD^2))^.5
  
  jointfullresult.m <- cbind(meancmaqp, combinedmu, combinedse)
  JointFullResults = rbind (JointFullResults, jointfullresult.m)

} #end new location loop
list(ensemble_weight = JointFullResults$meancmaqp, est_mean = JointFullResults$combinedmu2, est_se = JointFullResults$combinedse2)
}#end pred_spatial_ensemble function



#####################
#END ensemble weight functions
#####################


