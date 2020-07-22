# This function draws the fixed effects parameter beta marginally
#--------------------------------------------------------------------------
  # Output
# beta      : Gibbs draw of beta, a column vector (k x 1)
# Input
# z         : Latent response variable, matrix of size (m x n)
# x         : covariates including a column of ones, size (k x m x n)
# s         : random-effects covariates, size (l x m x n)
# alpha     : random-effects parameter has been integrated out
# w         : matrix, size (m x n)
# varphi2   : scalar quantity
# tau2      : parameter from normal-exponential mixture of Laplace dist.
# theta     : parameter from normal-exponential mixture of Laplace dist.
# invB0     : inverse of the prior covariance matrix
# invB0b0   : prior precision times prior mean
#--------------------------------------------------------------------------

sampleBeta <- function(z,x,s,w,varphi2,tau2,theta,invB0,invB0b0,k,m,n)
{
  sumvar  = matrix(0,nrow=k,ncol=k)
  summean = matrix(0,nrow=k,ncol=1)

  for(i in 1:n)
  {
    #D1 = tau2*(diag(w[,i]))
    inv_omega = x[,,i]%*%solve(varphi2*(t(s[,,i])%*%s[,,i]) + tau2*(diag(w[,i])))

    sumvar  = sumvar + inv_omega%*%t(x[,,i]) 
    summean = summean + inv_omega%*%(z[,i] - theta*w[,i])
  }

  Btilde = solve(invB0+sumvar)
  btilde  = Btilde%*%(invB0b0 + summean)

  beta = mvtnorm::rmvnorm(1,btilde,Btilde)
  return(beta)
}

# ALTERNATE without stacking is an option - check for matrix multiplication speeds.
