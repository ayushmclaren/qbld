#' @title Asymmetric Laplace distribution
#'
#'
#' @description Cumulative density, probability distribution function, quantile
#' function and random generation for the asymmetric Laplace distribution with
#' quantile \eqn{p}{p}, location parameter \code{mu} and scale parameter \code{sigma}.
#'
#' @details The asymmetric Laplace distribution (ALD), which has the following pdf:
#'
#'\deqn{f(x;\mu,\sigma,p) = \frac{p(1-p)}{\sigma} \exp\{-\frac{(x-\mu)}{\sigma}(p-I(x \le \mu))\}}{f(x;\mu,\sigma,p) = p(1-p)/\sigma exp{-(x-\mu)/\sigma (p-I(x \le \mu))}}
#'
#' If not specified, \eqn{p=0.5}, \eqn{mu = 0}, \eqn{sigma = 1}.
#'
#' @name aldmix
#' @aliases raldmix daldmix paldmix qaldmix
#' @param mu : location parameter
#' @param sigma : scale parameter
#' @param n : number of observations
#' @param x,q : vector of quantiles
#' @param p,prob       : probability at which to calculate quantile
#' @param lower.tail : logical; decides b/w \eqn{P(X<=p)} or \eqn{P(X>p)} for p/q
#'
#' @return
#' \itemize{
#' \item {\code{raldmix}} {returns a vector of random numbers from \code{AL(mu,sigma,p).}}
#' \item {\code{daldmix}} {returns returns density of \code{AL(mu,sigma,p)} at point x.}
#' \item {\code{paldmix}} {returns CDF prob of \code{AL(mu,sigma,p)} at quantile q.}
#' \item {\code{qaldmix}} {returns inverse CDF quantile of \code{AL(mu,sigma,p)} at prob.}
#' }
#'
#' @examples
#' raldmix(n = 10, mu = 5, sigma = 10, p = 0.5)
#' daldmix(c(4,5),mu = 0,sigma = 1,p = 0.5)
#' paldmix(c(1,4),mu = 0,sigma = 1,p = 0.5,lower.tail=TRUE)
#' qaldmix(0.5,mu = 0,sigma = 1,p = 0.5,lower.tail=TRUE)
#'
#' @references
#'
#' Keming Yu & Jin Zhang (2005) A Three-Parameter Asymmetric Laplace Distribution
#' and Its Extension, Communications in Statistics - Theory and Methods, 34:9-10,
#' 1867-1879, DOI: 10.1080/03610920500199018
#'
#' Kobayashi, Genya. (2011). Gibbs Sampling Methods for Bayesian Quantile Regression.
#' J Stat Comput Simul. 81. 1565. 10.1080/00949655.2010.496117.
#' 
#' @seealso \code{\link{rgig}} for random sampling from GIG distribution 


## usethis namespace: start
#' @useDynLib qbld, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @rdname aldmix
#' @export
raldmix <- function(n, mu, sigma, p) {
  .Call(`_qbld_raldmix`, n, mu, sigma, p)
}


#' @rdname aldmix
#' @export
daldmix = function(x,mu=0,sigma=1,p=0.5)
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

#' @rdname aldmix
#' @export
paldmix = function(q,mu=0,sigma=1,p=0.5,lower.tail=TRUE)
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

#' @rdname aldmix
#' @export
qaldmix = function(prob,mu=0,sigma=1,p=0.5,lower.tail=TRUE)
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


# Microbechmark test for speed/efficiency of rald_mix vs the existing implementation
# require(ald)
# require(ggplot2)
# mbm <- microbenchmark(my = rald_mix(1e5,0,1,0.5),old=rALD(1e5,0,1,0.5),times = 50)
# print(mbm)
# autoplot(mbm)
