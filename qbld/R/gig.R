#' @title Generalised Inverse Gaussian
#' @name gig
#' @aliases rgig dgig
#' 
#' @description Probability distribution function, Quantile Random generation 
#' for the Generalised Inverse Gaussian with three parameters a,b,p.
#' 
#' @details The Generalised Inverse Gaussian distrubtion(GIG), which has the following pdf
#' 
#' \deqn{f(x) = x^{\lambda-1}exp{-\frac{\omega}{2}(x + \frac{1}{x})}}{f(x) = x^(\lambda-1) exp{-\omega/2 (x + 1/x)}}
#' 
#' @param lambda,p : lambda parameter
#' @param a : chi parameter
#' @param b : psi parameter
#' @param n : number of observations
#' @param log_density : returns log density
#' @param x : vector of points
#' 
#' @return
#' \itemize{
#' \item {\code{rgig}} {returns a vector of random numbers from GIG(a,b,p).}
#' \item {\code{dgig}} {returns returns density of a GIG(a,b,p) at point x.}
#' }
#' 
#' @examples
#' rgig(1,0.5, 1, 2)
#' dgig(1, 1, 2, 0.5, FALSE)
#' 
#' @references
#' Devroye, L. Random variate generation for the generalized inverse Gaussian distribution. 
#' Stat Comput 24, 239–246 (2014).
#' 


#' @useDynLib qbld
#' @importFrom Rcpp sourceCpp
NULL

#'@rdname gig
#'@export
dgig <- function(x, a, b, p, log_density) {
  .Call(`_qbld_dgig`, x, a, b, p, log_density)
}

#' @useDynLib qbld
#' @importFrom Rcpp sourceCpp
NULL

#'@rdname gig
#'@export
rgig <- function(n, lambda, a, b) {
  .Call(`_qbld_rgig`, n, lambda, a, b)
}