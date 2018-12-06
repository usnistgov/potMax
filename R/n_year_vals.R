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
gumbelNYear <- function(x, N, ...) {
  UseMethod('gumbelNYear')
}

#' @export
gumbelNYear.gumbel_pot_fit <- function(x, N) {

  gumbelNYear.default(x = x$par,
                      thresh = x$thresh,
                      N = N)
}

#' @export
gumbelNYear.default <- function(x, thresh, N) {

  mu <- x[1]
  sigma <- x[2]

  N_year_val <- sigma*log(N) + mu

  value <- list(par = x,
                thresh = thresh,
                N = N,
                N_year_val = N_year_val)
  class(value) <- 'gumbel_N_year_val'
  value
}

#' @export
gumbelNYearUncert <- function(x, N, ...) {
  UseMethod('gumbelNYearUncert')
}

#' @export
gumbelNYearUncert.gumbel_pot_fit <- function(x, N, n_boot) {

  gumbelNYearUncert.default(x = x$par,
                            cov_mat = -solve(x$lhessian),
                            thresh = x$thresh,
                            N = N,
                            n_boot = n_boot)
}

#' @export
gumbelNYearUncert.default <- function(x, cov_mat,
                                      thresh, N,
                                      n_boot) {

  mu <- x[1]
  sigma <- x[2]

  bootstrap_samples <- MASS::mvrnorm(n = n_boot, mu = c(mu, log(sigma)), Sigma = cov_mat)
  bootstrap_samples[, 2] <- exp(bootstrap_samples[, 2])

  N_year_vals <- bootstrap_samples[, 2]*log(N) + bootstrap_samples[, 1]

  value <- list(par = c(mu, sigma),
                cov_mat = cov_mat,
                thresh = thresh,
                N = N,
                boot_samps = N_year_vals)
  class(value) <- 'gumbel_N_year_val_uncert'
  value
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
  UseMethod('fullNYear')
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
