#' @title Generalised Inverse Gaussian
#'
#' @aliases rgig dgig
#'
#' @description Probability distribution function, Quantile Random generation
#' for the Generalised Inverse Gaussian with three-parameters a,b,p.
#'
#' @details The Generalised Inverse Gaussian distrubtion(GIG), which has the following pdf:
#'
#' \deqn{f(x) = x^{\lambda -1}exp{-\frac{\omega}{2}(x + \frac{1}{x})}}{f(x) = x^(\lambda -1) exp{-\omega/2 (x + 1/x)}}
#'
#' @name GIGpack
#' @param p,P : lambda parameter
#' @param a : chi parameter
#' @param b : psi parameter
#' @param n : number of observations
#' @param log_density : returns log density
#' @param x : vector of points
#'
#' @return
#' \itemize{
#' \item {\code{rgig}} {returns a vector of random numbers from GIG(a,b,p).}
#' \item {\code{dgig}} {returns returns density of AGIG(a,b,p)at point x.}
#' }
#'
#' @examples
#' rgig (0.5, 1, 2, 1)
#' dgig(x=rep(1,2), 1, 2, 0.5,FALSE)
#'
#' @references
#'Devroye, L. Random variate generation for the generalized inverse Gaussian distribution.
#'Stat Comput 24, 239–246 (2014). https://doi.org/10.1007/s11222-012-9367-z
#'

#' @useDynLib GIGpack
#' @importFrom Rcpp sourceCpp
NULL
#' @rdname GIGpack
#' @export

dgig <- function(x, a, b, p, log_density) {
  .Call('_GIGpack_dgig', PACKAGE = 'GIGpack', x, a, b, p, log_density)
}

#' @useDynLib GIGpack
#' @importFrom Rcpp sourceCpp
NULL

#' @rdname GIGpack
#' @export

rgig <- function(P, a, b, n) {
  .Call('_GIGpack_rgig', PACKAGE = 'GIGpack', P, a, b, n)
}
