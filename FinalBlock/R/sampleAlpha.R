# This function samples the random effects parameter alpha marginally
#--------------------------------------------------------------------------
# Output
# alpha     : Gibbs draw of alpha, matrix of size (l x n)
# Input
# z         : Latent response variable, matrix of size (m x n)
# x         : covariates including a column of ones, size (k x m x n)
# s         : random-effects covariates, size (l x m x n)
# beta      : fixed-effects parameter, size (k x 1)
# w         : matrix, size (m x n)
# tau2      : parameter from normal-exponential mixture of Laplace dist.
# theta     : parameter from normal-exponential mixture of Laplace dist.
# varphi2   : scalar quantity
#--------------------------------------------------------------------------

sampleAlpha <- function(z,x,s,beta,w,tau2,theta,varphi2)
{
  l = dim(s)[1]
  m = dim(s)[2]
  n = dim(s)[3]

  alpha = matrix(0,nrow=l,ncol=n)
  invDvarphi2 = solve(varphi2*diag(l))

  for(i in 1:n)
  {
    invD1 = solve(tau2*diag(w[,i])) # do a 1/ instead of solve

    var = s[,,i]%*%invD1%*%t(s[,,i])
    mean = s[,,i]%*%invD1%*%(z[,i] - t(x[,,i])%*%beta - theta*w[,i]) #change names

    Atilde = solve(var + invDvarphi2)
    atilde = Atilde%*%(mean)
    alpha[,i] = mvtnorm::rmvnorm(1,mean=atilde,sigma=Atilde)
  }
  return(alpha)
}
