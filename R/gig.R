#' @title Generalised Inverse Gaussian
#' @name gig
#' @aliases rgig dgig
#' 
#' @description Probability distribution function, random generation 
#' for the Generalised Inverse Gaussian with three parameters \eqn{a(chi)}, \eqn{b(psi)}, \eqn{p}.
#' 
#' @details The Generalised Inverse Gaussian distrubtion(GIG), which has the following pdf
#' 
#' \deqn{f(x) = x^{\lambda-1}\exp\{-\frac{\omega}{2}(x + \frac{1}{x})\}}{f(x) = x^(\lambda-1) exp{-\omega/2 (x + 1/x)}}
#' 
#' @param lambda,p : lambda parameter
#' @param a : chi parameter. Must be nonnegative for positive lambda and positive else.
#' @param b : psi parameter. Must be nonnegative for negative lambda and positive else.
#' @param n : number of observations
#' @param log_density : logical; returns log density if TRUE
#' @param x : Argument of pdf
#' 
#' @return
#' \itemize{
#' \item {\code{rgig}} {returns a vector of random numbers from \code{GIG(a,b,p)}.}
#' \item {\code{dgig}} {returns returns density of a \code{GIG(a,b,p)} at point x.}
#' }
#' 
#' @examples
#' rgig(n = 1, lambda = 0.5, a = 1, b = 2)
#' dgig(x = 1, a = 1, b = 2, p = 0.5, log_density = FALSE)
#' 
#' @references
#' Devroye, L. Random variate generation for the generalized inverse Gaussian distribution. 
#' Stat Comput 24, 239–246 (2014).
#' 
#' Wolfgang Hörmann and Josef Leydold (2013). 
#' Generating generalized inverse Gaussian random variates, Statistics and Computing (to appear), 
#' DOI: 10.1007/s11222-013-9387-3
#'
#' J. S. Dagpunar (1989). An easily implemented generalised inverse Gaussian generator, 
#' Comm. Statist. B – Simulation Comput. 18, 703–710.
#' 
#' 
#' @seealso \code{\link{raldmix}} for random sampling from Asymmetric Laplace distribution 


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