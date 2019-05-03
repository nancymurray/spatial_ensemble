# Nancy Murray 
# Need corresponding sample data and functions from 
#https://github.com/nancymurray/spatial_ensemble

#Load the training data from your local drive
load("train.data")

#queue up functions from your local drive
source ("spatial_ensemble_fn.R")

#rename datasets to match default in the spatial_ensemble function
train.cmaq.dat <- other_clust_CMAQ 
train.aod.dat <- other_clust_AOD

#rename datasets to match default in the pred_spatial_ensemble function
test.cmaq.dat <- sing_clust_CMAQ 
test.aod.dat <- sing_clust_AOD

#estimation of weights
saved.wt <- spatial_ensemble(train.cmaq.dat=train.cmaq.dat, train.aod.dat=train.aod.dat, Dist.mat=Dist.mat,
         n.iter = 10000, n.burn = 2500,  rho.tune=0.4, tunevar=0.45, 
         tau2= 1, rho = 10, a = 0.001, b = 0.001, c = 5, d = 0.05)

#prediction of weights to other locations
pred_spatial_ensemble(obj=saved.wt, Dist.mat=Dist.mat, Dist.mat.new = Dist.mat.new, 
                      test.cmaq.dat=test.cmaq.dat, test.aod.dat=test.aod.dat)