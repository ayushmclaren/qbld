# This function samples the variance of random-effects, varphi2
#--------------------------------------------------------------------------
  # Output
# varphi2   : a scalar quantity
# Input
# alpha     : random-effects parameter, size (l x n)
# c1        : IG prior hyperparameter
# d1        : IG prior hyperparameter
#--------------------------------------------------------------------------

sampleVarphi2 <-function(alpha,c1,d1)
{
  l = nrow(alpha)
  n = ncol(alpha)
  scale  = matrix(0,nrow=n,ncol=1)

  for(i in 1:n)
  scale[i,1] = t(alpha[,i])%*%alpha[,i]


  tilde_c1 = (n*l+c1)/2;
  tilde_d1 = (sum(scale) + d1)/2;
  varphi2  = 1/rgamma(1,tilde_c1,tilde_d1)
    #1/rgig(tilde_c1,2*tilde_d1,0,1) #gamma(a,b) from gig(a,2b,0), inv_gamma = 1/gamma
}
#--------------------------------------------------------------------------
