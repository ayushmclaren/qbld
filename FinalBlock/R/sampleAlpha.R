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

sampleAlpha <- function(z,x,s,beta,w,tau2,theta,varphi2,l,m,n)
{
  alpha = matrix(0,nrow=l,ncol=n)
  invDvarphi2 = diag(l)/varphi2

  for(i in 1:n)
  {
    invD1 = s[,,i]%*%solve(tau2*diag(w[,i])) # do a 1/ instead of solve

    Atilde= solve(invD1%*%t(s[,,i]) + invDvarphi2)
    atilde = Atilde%*%(invD1%*%(z[,i] - t(x[,,i])%*%beta - theta*w[,i])) #change names
    alpha[,i] = mvntnorm::rmvnorm(1,atilde,Atilde)
  }
  return(alpha)
}
