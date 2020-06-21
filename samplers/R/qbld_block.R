##################################################################
### This is a skeleton, please DON'T RUN for now!
### DO NOT BUILD THIS PACKAGE, for internal purpose only!
##################################################################
### This essentially brings all the sampling steps together!
### The whole thing, with all the samplers will be migrated
### to Rcpp/Armadillo at a later stage for faster and efficient sampling!
##################################################################
### Convergence Diagnostics using mcmcse!
### as well as other outputs will be dealt with next!
##################################################################

#--------------------------------------------------------------------------
  ## GIBBS SAMPLING

tic
for nsim = 2:MCMC

#---------------------------------------------------------------------------------
  # The sequence of sampling is important and follows Algorithm 5 in Chib and Carlin.
#---------------------------------------------------------------------------------

  #---------- Sample beta,z marginally of alpha in a block --------------

  betap(:,nsim) = sampleBetaMarginal(z,x,s,w,varphi2(nsim-1,1),tau2,theta,invB0, invB0b0);

#---------- Sample z, marginally of alpha -----------------------------
  # Draws random numbers from trucnated MVN based on a series of conditional posteriors (Geweke, 1991).

zPrev = z;
z = sampleZMarginalCP(zPrev,y,x,s,betap(:,nsim),varphi2(nsim-1,1),w,tau2,theta,lowerLimits,upperLimits);

#------------- Sample alpha conditionally on beta,z -------------------
  alpha(:,:,nsim) = sampleAlpha(z,x,s,betap(:,nsim),w,tau2,theta,varphi2(nsim-1,1));

#---------- Sample w ---------------------------
  w = sampleW(z,x,s,betap(:,nsim),alpha(:,:,nsim),tau2,theta,lambdaIndex);

#---------- Sample varphi2 ---------------------
  varphi2(nsim,1) = samplevarphi2(alpha(:,:,nsim),c1,d1);

toc
