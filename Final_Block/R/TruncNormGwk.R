erf <- function(x, inverse = FALSE) {
  if (inverse) {
    ans <- qnorm((x+1)/2) / sqrt(2)
    ans[x <  -1] <- NA
    ans[x >  +1] <- NA
    ans[x == -1] <- -Inf
    ans[x == +1] <-  Inf
    ans
  } else {
    2 * pnorm(x * sqrt(2)) - 1
  }
}


erfc <- function(x, inverse = FALSE) {
  if (inverse) {
    ans <- qnorm(x/2, lower.tail = FALSE) / sqrt(2)
    ans[x <  0] <- NA
    ans[x >  2] <- NA
    ans[x == 0] <-  Inf
    ans[x == 2] <- -Inf
    ans
  } else {
    2 * pnorm(x * sqrt(2), lower.tail = FALSE)
  }
}


normcdf <- function(x)
{
  return(0.5*erfc(-x/sqrt(2)))
}

norminv <- function(x)
{
  return(-sqrt(2)*erfc(2*x,inverse = TRUE))
}

rtruncnorm_gwk <- function (z0,mu,sigma,LL,UL)
{
  
  m = length(z0)
  
  flag = 1
  if (nrow(sigma) != m) stop(paste("length(mean)=", m, ", but length(sigma)=", length(sigma), "!", sep=""))
  if(sum(z0[is.na(z0)==TRUE]) > 0) stop("NA's detected in z0.")
  if(sum(mu[is.na(mu)==TRUE]) > 0) stop("NA's detected in mu.")
  if(sum(UL[is.na(UL)==TRUE]) > 0) stop("NA's detected in UL.")
  if(sum(LL[is.na(LL)==TRUE]) > 0) stop("NA's detected in LL.")
  if (is.null(dim(sigma))) stop(paste("length(mean)=", m, ", but dim(sigma) is NULL!", sep=""))
  Q <- chol2inv(chol(sigma))       ## To check whether Sigma is PD
  
  if(any(!((LL[flag==1] <= z0[flag==1]) || (z0[flag==1] <= UL[flag==1])))) stop("z0 should be in LL to UL.")
  
  z = rep(0,m)
  newLL   = rep(0,m)
  newUL   = rep(0,m)
  condMu  = rep(0,m)
  condVar = rep(0,m)
  
  for (i in (1:m))
  {
    mu_i = mu[i]
    mu_less_i = mu[-i]
    
    sigma11 = sigma[i,i]
    sigma12 = t(as.matrix(sigma[i,-i]))
    sigma22 = sigma[-i,-i]
    
    inv_sigma22 = solve(sigma22)
    
    condVar[i] = sigma11 - sigma12 %*% inv_sigma22 %*% t(sigma12)
    
    z_less_i = z0[-i]
    
    condMu[i] = mu_i + sigma12 %*% inv_sigma22 %*% (z_less_i - mu_less_i)
    
    newLL[i] = (LL[i]- condMu[i])/sqrt(condVar[i])
    newUL[i] = (UL[i] - condMu[i])/sqrt(condVar[i])
    
    u = runif(1);
    z[i] =  condMu[i] + (sqrt(condVar[i]))*(norminv( u*(normcdf(newUL[i])
                         - normcdf(newLL[i])) + normcdf(newLL[i])))
    z0[i] = z[i]
  }
  
  return(z)
}