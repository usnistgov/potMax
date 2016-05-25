#' @title gumbelNYear
#'
#' @description gumbelNYear
#'
#' @details Calculates the N-year return value from a fitted POT model with the
#'   Gumbel tail length parameter, i.e., 0
#'
#' @param x
#'
#' @param N The N in N-year return value
#'
#' @export
#'
#' @useDynLib potMax
#'
#' @importFrom Rcpp evalCpp
#'
gumbelNYear <- function (x, N, ...) {
  UseMethod('gumbelNYear')
}

#' @export
gumbelNYear.gumbel_pot_fit <- function (x, N) {

  gumbelNYear.default(x = x$par,
                      N = N)
}

#' @export
gumbelNYear.default <- function (x, N) {

  x[2]*log(N) + x[1]
}

#' @title fullNYear
#'
#' @description fullNYear
#'
#' @details Calculates the N-year return value from a fitted POT model with
#'   a non-zero tail length parameter
#'
#' @param x
#'
#' @param N The N in N-year return value
#'
#' @export
#'
#' @useDynLib potMax
#'
#' @importFrom Rcpp evalCpp
#'
fullNYear <- function (x, N, ...) {
  UseMethod('gumbelNYear')
}

#' @export
fullNYear.full_pot_fit <- function (x, N) {

  fullNYear.default(x = x$par,
                    N = N)
}

#' @export
fullNYear.default <- function (x, N) {

  (((1/N)^(-x[3])) - 1)*(x[2]/x[3]) + x[1]
}