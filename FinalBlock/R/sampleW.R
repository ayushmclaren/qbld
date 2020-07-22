# This function samples the latent weight W
#--------------------------------------------------------------------------
  # Output
# w        : matrix of size (m x n) drawn from GIG distribution
# Input
# z         : Latent response variable, matrix of size (m x n)
# x         : covariates including a column of ones, size (k x m x n)
# s         : random-effects covariates, size (k x m x n)
# beta      : fixed-effects parameter, size (m x 1)
# alpha     : random-effects parameter, size (l x n)
# tau2      : paramter from normal-exponential mixture of Laplace dist.
# theta     : paramter from normal-exponential mixture of Laplace dist.
# lambda    : GIG parameter
#--------------------------------------------------------------------------

#sourceCpp("./src/rgig.cpp")

sampleW <-function(z,x,s,beta,alpha,tau2,theta,lambda=0.5,k,m,n)
{

  tilde_eta = (theta^2)/(tau2) + 2
  #tilde_lambda = matrix(0,nrow=m,ncol=n)
  w =  matrix(0,nrow=m,ncol=n)

  for(i in 1:n)
  {
    for(j in 1:m)
    {
      tilde_lambda = ((z[j,i] - t(x[,j,i])%*%beta - t(s[,j,i])%*%alpha[,i])^2)/tau2
      #print(tilde_lambda)
      
      if(tilde_lambda < 1e-8||is.nan(tilde_lambda)) tilde_lambda = 1e-8

      w[j,i] = rgig(lambda,tilde_eta,tilde_lambda,1)
    }
  }
  return(w)
}
