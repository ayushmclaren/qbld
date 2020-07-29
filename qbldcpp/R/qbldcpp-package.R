## usethis namespace: start
#' @useDynLib qbldcpp, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @rdname qbldcpp
#' @export
qbldf <- function(nsim, p, y, datax, datas, x_intercept, s_intercept, b0, B0, c1, d1, burnin) {
  .Call(`_qbldcpp_qbldf`, nsim, p, y, datax, datas, x_intercept, s_intercept, b0, B0, c1, d1, burnin)
}

## usethis namespace: start
#' @useDynLib qbldcpp, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @rdname qbldcpp
#' @export
datagen <- function(n, m, p) {
  .Call(`_qbldcpp_datagen`, n, m, p)
}

## usethis namespace: start
#' @useDynLib qbldcpp, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @rdname qbldcpp
#' @export
raldmix <- function(n, mu, sigma, p) {
  .Call(`_qbldcpp_raldmix`, n, mu, sigma, p)
}

## usethis namespace: start
#' @useDynLib qbldcpp, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @rdname qbldcpp
#' @export
rgig <- function(n, lambda, a, b) {
  .Call(`_qbldcpp_rgig`, n, lambda, a, b)
}

## usethis namespace: start
#' @useDynLib qbldcpp, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @rdname qbldcpp
#' @export
qbldunblock <- function(nsim, p, y, datax, datas, x_intercept, s_intercept, b0, B0, c1, d1, burnin) {
  .Call(`_qbldcpp_qbldunblock`, nsim, p, y, datax, datas, x_intercept, s_intercept, b0, B0, c1, d1, burnin)
}

