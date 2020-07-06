# function [z] = sampleZMarginalCP(zPrev,y,x,s,betap,varphi2,w,tau2,theta,lowerLimits,upperLimits)

# This function samples the latent variable Z marginally
#--------------------------------------------------------------------------
  # Output
# z        : matrix of size (m x n) latent variables
# Input
# zPrev     : Previous draw in the MCMC
# y         : observed response variable, matrix of size (m x n)
# x         : covariates including a column of ones, size (k x m x n)
# s         : random-effects covariates, size (l x m x n)
# beta      : fixed-effects parameter, size (k x 1)
# varphi2   : scalar, variance of random effects
# w         : matrix, size (m x n)
# tau2      : parameter from normal-exponential mixture of Laplace dist.
# theta     : parameter from normal-exponential mixture of Laplace dist.
# lowerLimits: lower truncation points, size m x n
# upperLimits: upper truncation points, size m x n
#--------------------------------------------------------------------------
#source("./R/TruncNormGwk.R")

sampleZ <- function(zprev,y,x,beta,s,theta,w,varphi2,tau2,LL,UL)
{

  m = nrow(y)
  n = ncol(y)

  z = matrix(0,nrow=m,ncol=n)

  for(i in 1:n)
  {
    meani = t(x[,,i])%*%beta + theta*w[,i]
    VarZi = (varphi2*(t(s[,,i])%*%s[,,i]) + tau2*(diag(w[,i],nrow=m)))

    z[,i] = rtruncnorm_gwk(zprev[,i],meani,VarZi,LL[,i],UL[,i])
  }
  return(z)
}
