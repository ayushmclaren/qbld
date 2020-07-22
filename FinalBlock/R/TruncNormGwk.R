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


limits <- function(y,which)
{
  if(which == 0)  #lowerlimits
  {
    if(y == 0) return(-Inf)
    else 
      return(0)
  }
  
  else if(which==1) #upperlimits
  {
    if(y == 0) return(0)
    else 
      return(Inf)
  }
}



rtruncnorm_gwk <- function (z0,mu,sigma,y,m,idx)
{
  
  #m = length(z0)
  flag = 1
  if (nrow(sigma) != m) stop(paste("length(mean)=", m, ", but length(sigma)=", length(sigma), "!", sep=""))
  if(sum(z0[is.na(z0)==TRUE]) > 0) stop("NA's detected in z0.")
  if(sum(mu[is.na(mu)==TRUE]) > 0) stop("NA's detected in mu.")
  #if(sum(UL[is.na(UL)==TRUE]) > 0) stop("NA's detected in UL.")
  #if(sum(LL[is.na(LL)==TRUE]) > 0) stop("NA's detected in LL.")
  if (is.null(dim(sigma))) stop(paste("length(mean)=", m, ", but dim(sigma) is NULL!", sep=""))
  Q <- chol2inv(chol(sigma))       ## To check whether Sigma is PD
  
  #if(any(!((LL[flag==1] <= z0[flag==1]) || (z0[flag==1] <= UL[flag==1])))) stop("z0 should be in LL to UL.")
  
  for (i in (1:m))
  {
    mu_i = mu[i]
    mu_less_i = mu[-i]
    
    sigma11 = sigma[i,i]
    sigma12 = t(as.matrix(sigma[i,-i]))
    inv_sigma22 = solve(sigma[-i,-i])
    
    condVar = sigma11 - sigma12 %*% inv_sigma22 %*% t(sigma12)
    
    z_less_i = z0[-i]
    
    condMu = mu_i + sigma12 %*% inv_sigma22 %*% (z_less_i - mu_less_i)
    
    newLL = (limits(y[i,idx],0)- condMu)/sqrt(condVar)
    newUL = (limits(y[i,idx],1) - condMu)/sqrt(condVar)
  
    z0[i] =  condMu + (sqrt(condVar))*(norminv(runif(1)*(normcdf(newUL) - normcdf(newLL)) + normcdf(newLL)))
  }
  return(z0)
}