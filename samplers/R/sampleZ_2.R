
# This function samples the latent variable Z marginally
#--------------------------------------------------------------------------
# Output
# z        : matrix of size (m x n) latent variables
# Input
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

sampleZ_2 <- function(y,x,s,betap,alpha,varphi2,w,tau2,theta,lowerLimits,upperLimits)
{

  m = nrow(y)
  n = ncol(y)

  z = matrix(0,nrow=m,ncol=n)

  for(i in 1:n)
  {
    for(j in 1:m)
    {
      mean_comp = t(x[,j,i])%*%beta + t(s[,j,i])%*%alpha[,i] + theta*w[j,i]
      var_comp  = tau2*w[j,i]

      if(y[j,i] == 0)
      z[j,i] = truncnorm::rtruncnorm(1,a=-Inf,b=0,mean = mean_comp, sd = sqrt(var_comp))

      if(y[j,i] == 1)
        z[j,i] = truncnorm::rtruncnorm(1,a=0,b=Inf,mean = mean_comp, sd = sqrt(var_comp))
    }
  }
  return(z)
}
