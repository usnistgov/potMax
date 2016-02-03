#' @title gumbelMaxDist
#'
#' @description gumbelMaxDist
#'
#' @details Generates the empirical distribution of the maximum value from a
#'   Gumbel (zero tail length) POT model over a specified length of time
#'
#' @param mu Location parameter
#'
#' @param sigma Scale parameter
#'
#' @param thresh The threshold
#'
#' @param lt_gen Lengt of each generated series in seconds
#'
#' @export
#'
#' @useDynLib potMax
#'
#' @importFrom Rcpp evalCpp
#'
gumbelMaxDist <- function (x, lt_gen, n_mc, ...) {
  UseMethod('gumbelMaxDist')
}

#' @export
gumbelMaxDist.gumbel_pot_fit <- function (x, lt_gen, n_mc) {

  gumbelMaxDist.default(x = x$par,
                        thresh = x$thresh,
                        lt_gen = lt_gen,
                        n_mc = n_mc)
}

#' @export
gumbelMaxDist.default <- function (x, thresh,
                                   lt_gen, n_mc) {

  mu <- x[1]
  sigma <- x[2]
  lt_gen <- lt_gen[1] # only scalars are allowed

  if (sigma <= 0) {
    stop('must have sigma > 0')
  }
  if (lt_gen <= 0) {
    stop('must have lt_gen > 0')
  }

  Lambda <- lt_gen*exp((-(thresh - mu))/sigma)
  const <- Lambda/lt_gen

  prob_zero_obs <- dpois(0, Lambda)
  if (prob_zero_obs > 0.9) {

    warning(paste0('The probability of zero observations is ', prob_zero_obs, '\n'))
  }

  value <- list(par = c(mu, sigma),
                thresh = thresh,
                lt_gen = lt_gen,
                max_dist = gumbelMaxDistCpp(mu, sigma,
                                            Lambda, const, n_mc))
  class(value) <- 'gumbel_max_dist'
  value
}

#' @title gumbelMaxDistUncert
#'
#' @description gumbelMaxDistUncert
#'
#' @details Generates bootstrap samples of the empirical distribution of the
#'   maximum value from a Gumbel (zero tail length) POT model over a specified
#'   #'   length of time
#'
#' @param mu Location parameter
#'
#' @param sigma Scale parameter
#'
#' @param cov_mat The estimated covariance matrix for the estimation of mu and
#'   sigma
#'
#' @param thresh The threshold
#'
#' @param lt_gen Lengt of each generated series in seconds that the maximum is
#'   take over
#'
#' @param n_mc The number of Monte Carlo samples from which the distribution of
#'   the maximum is built
#'
#' @param n_boot The number of bootstrap samples to base uncertainties on
#'
#' @export
#'
#' @useDynLib potMax
#'
#' @importFrom Rcpp evalCpp
#'
gumbelMaxDistUncert <- function (x, lt_gen,
                                 n_mc,
                                 n_boot,
                                 ...) {
  UseMethod('gumbelMaxDistUncert')
}

#' @export
gumbelMaxDistUncert.gumbel_pot_fit <- function (x, lt_gen,
                                                n_mc,
                                                n_boot) {

  gumbelMaxDistUncert.default(x = x$par,
                              cov_mat = -solve(x$hessian),
                              thresh = x$thresh,
                              lt_gen = lt_gen,
                              n_mc = n_mc,
                              n_boot = n_boot)
}

#' @export
gumbelMaxDistUncert.default <- function (x, cov_mat,
                                         thresh, lt_gen,
                                         n_mc,
                                         n_boot) {

  mu <- x[1]
  sigma <- x[2]
  lt_gen <- lt_gen[1]

  if (sigma <= 0) {
    stop('must have sigma > 0')
  }
  if (lt_gen <= 0) {
    stop('must have lt_gen > 0')
  }

  bootstrap_samples <- MASS::mvrnorm(n = n_boot, mu = c(mu, sigma), Sigma = cov_mat)
  n_neg <- sum(bootstrap_samples[, 2] <= 0)
  if (n_neg > 0){

    bootstrap_samples <- bootstrap_samples[bootstrap_samples[, 2] > 0, ]
    warning(paste0('Only using ', n_boot - n_neg, ' bootstrap samples instead of ', n_boot))
    n_boot <- n_boot - n_neg
  }

  Lambda <- lt_gen*exp((-(thresh - bootstrap_samples[, 1]))/
                         bootstrap_samples[, 2])
  const <- Lambda/lt_gen

  prob_zero_obs <- dpois(0, Lambda)
  if (sum(prob_zero_obs > 0.9) > 0) {

    warning('The probability of zero observations for some of the bootstrap samples of the parameters is greater than 0.9')
  }

  value <- list(par = c(mu, sigma),
                cov_mat = cov_mat,
                thresh = thresh,
                lt_gen = lt_gen,
                boot_samps = gumbelMaxDistUncertCpp(bootstrap_samples[, 1],
                                                    bootstrap_samples[, 2],
                                                    Lambda, const,
                                                    n_mc, n_boot))
  class(value) <- 'gumbel_max_dist_uncert'
  value
}

#' @title fullMaxDist
#'
#' @description fullMaxDist
#'
#' @details Generates the empirical distribution of the maximum value from a
#'   full (non-zero tail length) POT model over a specified length of time
#'
#' @param mu Location parameter
#'
#' @param sigma Scale parameter
#'
#' @param k Tail length parameter
#'
#' @param thresh The threshold
#'
#' @param lt_gen Lengt of each generated series in seconds
#'
#' @export
#'
#' @useDynLib potMax
#'
#' @importFrom Rcpp evalCpp
#'
fullMaxDist <- function (x, lt_gen, n_mc, ...) {
  UseMethod('fullMaxDist')
}

#' @export
fullMaxDist.full_pot_fit <- function (x, lt_gen, n_mc) {

  fullMaxDist.default(x = x$par,
                      thresh = x$thresh,
                      lt_gen = lt_gen,
                      n_mc = n_mc)
}

#' @export
fullMaxDist.default <- function (x, thresh,
                                 lt_gen, n_mc) {

  mu <- x[1]
  sigma <- x[2]
  k <- x[3]
  lt_gen <- lt_gen[1] # only scalars are allowed

  if (sigma <= 0) {
    stop('must have sigma > 0')
  }
  if (lt_gen <= 0) {
    stop('must have lt_gen > 0')
  }

  Lambda <- 1 + k*(thresh - mu)/sigma
  Lambda <- Lambda^(-1/k)
  Lambda <- lt_gen*Lambda
  const <- Lambda/lt_gen

  prob_zero_obs <- dpois(0, Lambda)
  if (prob_zero_obs > 0.9) {

    warning(paste0('The probability of zero observations is ', prob_zero_obs, '\n'))
  }

  value <- list(par = c(mu, sigma, k),
                thresh = thresh,
                lt_gen = lt_gen,
                max_dist = fullMaxDistCpp(mu, sigma, k,
                                          Lambda, const, n_mc))
  class(value) <- 'full_max_dist'
  value
}

#' @title fullMaxDistUncert
#'
#' @description fullMaxDistUncert
#'
#' @details Generates bootstrap samples of the empirical distribution of the
#'   maximum value from the full (non-zero tail length) POT model over a
#'   specified length of time
#'
#' @param mu Location parameter
#'
#' @param sigma Scale parameter
#'
#' @param k Tail length parameter
#'
#' @param cov_mat The estimated covariance matrix for the estimation of mu and
#'   sigma
#'
#' @param thresh The threshold
#'
#' @param lt_gen Lengt of each generated series in seconds that the maximum is
#'   taken over
#'
#' @param n_mc The number of Monte Carlo samples from which the distribution of
#'   the maximum is built
#'
#' @param n_boot The number of bootstrap samples to base uncertainties on
#'
#' @export
#'
#' @useDynLib potMax
#'
#' @importFrom Rcpp evalCpp
#'
fullMaxDistUncert <- function (x, lt_gen,
                               n_mc,
                               n_boot,
                               ...) {
  UseMethod('fullMaxDistUncert')
}

#' @export
fullMaxDistUncert.full_pot_fit <- function (x, lt_gen,
                                            n_mc,
                                            n_boot) {

  fullMaxDistUncert.default(x = x$par,
                            cov_mat = -solve(x$hessian),
                            thresh = x$thresh,
                            lt_gen = lt_gen,
                            n_mc = n_mc,
                            n_boot = n_boot)
}

#' @export
fullMaxDistUncert.default <- function (x,
                                       cov_mat, thresh, lt_gen,
                                       n_mc,
                                       n_boot) {

  mu <- x[1]
  sigma <- x[2]
  k <- x[3]
  lt_gen <- lt_gen[1]

  if (sigma <= 0) {
    stop('must have sigma > 0')
  }
  if (lt_gen <= 0) {
    stop('must have lt_gen > 0')
  }

  bootstrap_samples <- MASS::mvrnorm(n = n_boot, mu = c(mu, sigma, k),
                                     Sigma = cov_mat)
  n_neg <- sum(bootstrap_samples[, 2] <= 0)
  if (n_neg > 0){
    bootstrap_samples <- bootstrap_samples[bootstrap_samples[, 2] > 0, ]
  }

  Lambda <- 1 + bootstrap_samples[, 3]*
    (thresh - bootstrap_samples[, 1])/bootstrap_samples[, 2]
  Lambda <- Lambda^(-1/bootstrap_samples[, 3])
  Lambda <- lt_gen*Lambda
  const <- Lambda/lt_gen

  bootstrap_samples <- bootstrap_samples[!is.na(Lambda), ]
  const <- const[!is.na(Lambda)]
  Lambda <- Lambda[!is.na(Lambda)]
  new_n_boot <- length(Lambda)
  if (new_n_boot < n_boot) {
    warning(paste0('Only using ', new_n_boot,
                   ' bootstrap samples instead of ',
                   n_boot))
  }
  prob_zero_obs <- dpois(0, Lambda)
  if (sum(prob_zero_obs > 0.9) > 0) {

    warning('The probability of zero observations for some of the bootstrap samples of the parameters is greater than 0.9')
  }

  value <- list(par = c(mu, sigma, k),
                cov_mat = cov_mat,
                thresh = thresh,
                lt_gen = lt_gen,
                boot_samps = fullMaxDistUncertCpp(bootstrap_samples[, 1],
                                                  bootstrap_samples[, 2],
                                                  bootstrap_samples[, 3],
                                                  Lambda, const,
                                                  n_mc, new_n_boot))
  class(value) <- 'full_max_dist_uncert'
  value
}
