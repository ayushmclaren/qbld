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


sampleAlphafast <- function(z,x,s,beta,w,tau,theta,varphi)
{
  l = dim(s)[1]
  m = dim(s)[2]
  n = dim(s)[3]
  
  out = matrix(0,nrow=l,ncol=n)
  # D = varphi2*diag(l)
  # phi = invD1 (sqrt) * S'
  # n = m, p = l
  # alpha = invD1(sqrt)*(zi - xi*beta - theta*wi)
  
  for(i in 1:n)
  {
    invD1 = solve(tau*diag(sqrt(w[,i])))
    phi = invD1%*%t(s[,,i])
    u = rnorm(l,0,varphi)
    delta = rnorm(m,0,1) 
    
    v = phi%*%u + delta
    alph = invD1%*%((z[,i] - t(x[,,i])%*%beta - theta*w[,i]))
    w_i = solve((varphi^2)*phi%*%t(phi) + diag(m), alph - v)
    out[,i] = u + (varphi^2)*t(phi)%*%w_i
  }
  return(out)
}