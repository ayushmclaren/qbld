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


sampleAlpha_fast <- function(z,x,s,beta,w,tau2,theta,varphi2,l,m,n)
{
  
  out = matrix(0,nrow=l,ncol=n)
  # D = varphi2*diag(l)
  # phi = invD1 (sqrt) * S'
  # n = m, p = l
  # alpha = invD1(sqrt)*(zi - xi*beta - theta*wi)
  
  for(i in 1:n)
  {
    invD1 = diag(1/sqrt(w[,i]*tau2))
    phi = invD1%*%t(s[,,i])
    u = rnorm(l,0,sqrt(varphi2))
    delta = rnorm(m,0,1) 
    
    v = phi%*%u + delta
    alph = invD1%*%((z[,i] - t(x[,,i])%*%beta - theta*w[,i]))
    w_i = solve((varphi2)*phi%*%t(phi) + diag(m), alph - v)
    out[,i] = u + (varphi2)*t(phi)%*%w_i
  }
  return(out)
}