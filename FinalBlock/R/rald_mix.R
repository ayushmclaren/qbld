
rald_mix <- function(n,mu=0,sigma=1,p=0.5,plot=FALSE)
{
  if(length(n) != 1) stop("Provide integer sample size n.")
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
