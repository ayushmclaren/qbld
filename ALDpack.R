#' These functions are the generic r/d/p/q implementations of
#' AL distribution using normal-exponential mixture representation.
#' --------------------------------------------------------------------------

#' This function generates a vector of random numbers from an asymmetric 
#' Laplace distribution with quantile p.
#' --------------------------------------------------------------------------
#'  Output    
#' r    : a vector (n x 1) of random numbers from an AL(mu,sigma,p)
#' Input
#' mu    : location parameter
#' sigma : scale parameter
#' p     : quantile (skewness parameter)
#' n     : number of observations
#' plot  : creates basic plot
#'--------------------------------------------------------------------------

rald_mix <- function(n,mu=0,sigma=1,p=0.5,plot=FALSE)  
{
  if(length(n) == 0) stop("Provide sample size n.")
  if(n <= 0 || n%%1 !=0) stop("Sample size must be a positive integer value.")
  if(sigma <= 0) stop("sigma (scale parameter) must be a positive number.")
  if(p >= 1 || p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu (location parameter) must be a finite real number.")
  
  theta <- (1-2*p)/(p*(1-p))
  tau <- sqrt(2/(p*(1-p)))
  
  z <- rexp(n,rate=1)
  u <- rnorm(n,mean=0,sd=1) 
  
  r <- theta*z + tau* sqrt(z) * u
  r <- r + rep(mu,n)
  
  if(plot == TRUE) #creates a sample plot for the ALD
  plot(density(rald_mix(1e5,mu,sigma,p)),main=paste0("ALD(", mu,",",sigma,",",p,")"),col="blue")
  
  return(r)
}

#' This function returns density of AL(mu,sigma,p)
#' --------------------------------------------------------------------------
#'  Output    
#' d   : density of AL(mu,sigma,p) at x
#' Input (additional)
#' x     : point at which to calculate 
#'--------------------------------------------------------------------------

dald_mix = function(x,mu=0,sigma=1,p=0.5)
{
  if(length(x) == 0) stop("Provide points x.")
  if(sum(x[is.na(x)==TRUE]) > 0) stop("NA's detected in x.")
  if(sigma <= 0) stop("sigma (scale parameter) must be a positive number.")
  if(p >= 1 || p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) == Inf) stop("mu (location parameter) must be a finite real number.")
  
  d = ifelse(test= (x < mu), yes = (p*(1-p)/sigma)*exp((1-p)*(x-mu)/sigma),
             no = (p*(1-p)/sigma)*exp(-(p)*(x-mu)/sigma))
  return(d)
}

#' This function returns CDF prob of AL(mu,sigma,p)
#' --------------------------------------------------------------------------
#'  Output    
#' pald   : probability
#' Input (additional)
#' q          : quantile
#' lower.tail : decides b/w P(X<=p) or P(X>p) 
#'--------------------------------------------------------------------------

pald_mix = function(q,mu=0,sigma=1,p=0.5,lower.tail=TRUE)
{
  if(length(q) == 0) stop("Provide quantile q.")
  if(sigma <= 0) stop("sigma (scale parameter) must be a positive number.")
  if(p >= 1 || p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu (location parameter) must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be specified.")
  
  pald = ifelse(test= (q < mu), yes = p*exp((1-p)*(q-mu)/sigma),
                no = 1-(1-p)*exp(-(p)*(q-mu)/sigma))
  
  ifelse(test=lower.tail == TRUE,yes=return(pald),no=return(1-pald))
  
} 

#' This function returns inverse CDF quantile of AL(mu,sigma,p)
#' --------------------------------------------------------------------------
#' Output    
#' qald   : quantile
#' Input (additional)
#' prob       : probability
#' lower.tail : decides b/w P(X<=p) or P(X>p) 
#'--------------------------------------------------------------------------

qald_mix = function(prob,mu=0,sigma=1,p=0.5,lower.tail=TRUE)
{
  if(length(prob) == 0) stop("Provide prob.")
  if(sigma <= 0) stop("sigma (scale parameter) must be a positive number.")
  if(p >= 1 || p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu (location parameter) must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be specified.")
  if(sum(prob > 1 | prob < 0) > 0) stop("All elements of prob must be real numbers in (0,1).")
  
  q = sapply(X=prob,FUN=function(prob,mu,sigma,p)
    
    {ifelse(test= (lower.tail == TRUE), yes = prob, no = (1-prob));
     ifelse(test= (prob<p), yes = mu+ ((sigma*log(prob/p))/(1-p)),
           no = mu - ((sigma*log((1-prob)/(1-p)))/p))}, mu=mu,sigma=sigma,p=p)
  
  return(q)
}


#' Microbechmark test for speed/efficiency of rald_mix vs the existing implementation
#' require(ald)
#' require(ggplot2)
#' mbm <- microbenchmark(my = rald_mix(1e5,0,1,0.5),old=rALD(1e5,0,1,0.5),times = 50)
#' print(mbm)
#' autoplot(mbm)
