## usethis namespace: start
#' @useDynLib qbild, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @rdname qbild
#' @export
qbldf <- function(nsim, p, y, datax, datas, x_intercept, s_intercept, b0, B0, c1, d1, m, n, k, l) {
  .Call(`_qbild_qbldf`, nsim, p, y, datax, datas, x_intercept, s_intercept, b0, B0, c1, d1, m, n, k, l)
}

## usethis namespace: start
#' @useDynLib qbild, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @rdname qbild
#' @export
datagen <- function(n, m, p) {
  .Call(`_qbild_datagen`, n, m, p)
}

## usethis namespace: start
#' @useDynLib qbild, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @rdname qbild
#' @export
raldmix <- function(n, mu, sigma, p) {
  .Call(`_qbild_raldmix`, n, mu, sigma, p)
}

## usethis namespace: start
#' @useDynLib qbild, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @rdname qbild
#' @export
rgig <- function(n, lambda, a, b) {
  .Call(`_qbild_rgig`, n, lambda, a, b)
}

## usethis namespace: start
#' @useDynLib qbild, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @rdname qbild
#' @export
qbldunblock <- function(nsim, p, y, datax, datas, x_intercept, s_intercept, b0, B0, c1, d1, m, n, k, l) {
  .Call(`_qbild_qbldunblock`, nsim, p, y, datax, datas, x_intercept, s_intercept, b0, B0, c1, d1, m, n, k, l)
}
