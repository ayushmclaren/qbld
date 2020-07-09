

#' @useDynLib FinalBlock, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
#' @rdname FinalBlock
#' @export

dgig <- function(x, a, b, p, log_density) {
  .Call(`_FinalBlock_dgig`, x, a, b, p, log_density)
}

#' @useDynLib FinalBlock, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
#' @rdname FinalBlock
#' @export

rgig <- function(P, a, b, n) {
  .Call(`_FinalBlock_rgig`, P, a, b, n)
}

