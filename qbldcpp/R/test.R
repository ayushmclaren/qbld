if(0)
{

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


sampleAlphafast <- function(z,x,s,beta,w,tau,theta,varphi)
{
  l = dim(s)[1]
  m = dim(s)[2]
  n = dim(s)[3]
  
  out = matrix(0,nrow=l,ncol=n)
  u = rnorm(l,0,varphi)
  delta = rnorm(m,0,1) 
  UIm = diag(m)
  # D = varphi2*diag(l)
  # phi = invD1 (sqrt) * S'
  # n = m, p = l
  # alpha = invD1(sqrt)*(zi - xi*beta - theta*wi)
  
  for(i in 1:n)
  {
    invD1 = solve(tau*diag(sqrt(w[,i])))
    phi = invD1%*%t(s[,,i])
    v = phi%*%u + delta
    alph = invD1%*%((z[,i] - t(x[,,i])%*%beta - theta*w[,i]))
    w_i = solve((varphi^2)*phi%*%t(phi) + UIm, alph - v)
    out[,i] = u + (varphi^2)*t(phi)%*%w_i
  }
  return(out)
}

}