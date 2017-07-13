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
gumbelMaxDist <- function(x, lt_gen, n_mc, progress_tf = TRUE, ...) {
  UseMethod('gumbelMaxDist')
}

#' @export
gumbelMaxDist.gumbel_pot_fit <- function(x, lt_gen, n_mc,
                                         progress_tf = TRUE) {

  gumbelMaxDist.default(x = x$par,
                        thresh = x$thresh,
                        lt_gen = lt_gen,
                        n_mc = n_mc,
                        progress_tf = progress_tf)
}

#' @export
gumbelMaxDist.default <- function(x, thresh,
                                  lt_gen, n_mc,
                                  progress_tf = TRUE) {

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

    warning(paste0('The probability of zero observations is ', prob_zero_obs, '\n'), immediate. = TRUE)
  }

  value <- list(par = c(mu, sigma),
                thresh = thresh,
                lt_gen = lt_gen,
                max_dist = gumbelMaxDistCpp(mu, sigma,
                                            Lambda, const, n_mc,
                                            progress_tf))
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
gumbelMaxDistUncert <- function(x, lt_gen,
                                n_mc,
                                n_boot,
                                progress_tf = TRUE,
                                ...) {
  UseMethod('gumbelMaxDistUncert')
}

#' @export
gumbelMaxDistUncert.gumbel_pot_fit <- function(x, lt_gen,
                                               n_mc,
                                               n_boot,
                                               progress_tf = TRUE) {

  gumbelMaxDistUncert.default(x = x$par,
                              cov_mat = -solve(x$hessian),
                              thresh = x$thresh,
                              lt_gen = lt_gen,
                              n_mc = n_mc,
                              n_boot = n_boot,
                              progress_tf = progress_tf)
}

#' @export
gumbelMaxDistUncert.default <- function(x, cov_mat,
                                        thresh, lt_gen,
                                        n_mc,
                                        n_boot,
                                        progress_tf = TRUE) {

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
  if (n_neg > 0) {

    bootstrap_samples <- bootstrap_samples[bootstrap_samples[, 2] > 0, ]
    warning(paste0('Removing ', n_neg, ' bootstrap samples because sigma^* <= 0'), immediate. = TRUE)
    n_boot <- n_boot - n_neg
  }

  Lambda <- lt_gen*exp((-(thresh - bootstrap_samples[, 1])) /
                         bootstrap_samples[, 2])
  Lambda_orig <- lt_gen*exp((-(thresh - mu))/sigma)
  const <- Lambda/lt_gen

  prob_zero_obs <- dpois(0, Lambda)
  zero_prob_high <- prob_zero_obs > 0.9
  if (sum(zero_prob_high) > 0) {

    bootstrap_samples <- bootstrap_samples[!zero_prob_high, ]
    Lambda <- Lambda[!zero_prob_high]
    const <- const[!zero_prob_high]
    warning(paste0('Removing ', sum(zero_prob_high), ' bootstrap samples because the probability of zero threshold exceedances is > 90%'), immediate. = TRUE)
    n_boot <- n_boot - sum(zero_prob_high)
  }

  big_Lambda <- Lambda > 10*Lambda_orig
  if (sum(big_Lambda) > 0) {

    bootstrap_samples <- bootstrap_samples[!big_Lambda, ]
    Lambda <- Lambda[!big_Lambda]
    const <- const[!big_Lambda]
    warning(paste0('Removing ', sum(big_Lambda), ' bootstrap samples because Lambda* > 10Lambda_orig'), immediate. = TRUE)
    n_boot <- n_boot - sum(big_Lambda)
  }

  value <- list(par = c(mu, sigma),
                cov_mat = cov_mat,
                thresh = thresh,
                lt_gen = lt_gen,
                boot_samps = gumbelMaxDistUncertCpp(bootstrap_samples[, 1],
                                                    bootstrap_samples[, 2],
                                                    Lambda, const,
                                                    n_mc, n_boot, progress_tf))
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
fullMaxDist <- function(x, lt_gen, n_mc,
                        progress_tf = TRUE, ...) {
  UseMethod('fullMaxDist')
}

#' @export
fullMaxDist.full_pot_fit <- function(x, lt_gen, n_mc,
                                     progress_tf = TRUE) {

  fullMaxDist.default(x = x$par,
                      thresh = x$thresh,
                      lt_gen = lt_gen,
                      n_mc = n_mc,
                      progress_tf = progress_tf)
}

#' @export
fullMaxDist.default <- function(x, thresh,
                                lt_gen, n_mc,
                                progress_tf = TRUE) {

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

    warning(paste0('The probability of zero observations is ', prob_zero_obs, '\n'), immediate. = TRUE)
  }

  value <- list(par = c(mu, sigma, k),
                thresh = thresh,
                lt_gen = lt_gen,
                max_dist = fullMaxDistCpp(mu, sigma, k,
                                          Lambda, const, n_mc,
                                          progress_tf))
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
fullMaxDistUncert <- function(x, lt_gen,
                              n_mc,
                              n_boot,
                              progress_tf = TRUE,
                              ...) {
  UseMethod('fullMaxDistUncert')
}

#' @export
fullMaxDistUncert.full_pot_fit <- function(x, lt_gen,
                                           n_mc,
                                           n_boot,
                                           progress_tf = TRUE) {

  fullMaxDistUncert.default(x = x$par,
                            cov_mat = -solve(x$hessian),
                            thresh = x$thresh,
                            lt_gen = lt_gen,
                            n_mc = n_mc,
                            n_boot = n_boot,
                            progress_tf = progress_tf)
}

#' @export
fullMaxDistUncert.default <- function(x,
                                      cov_mat, thresh, lt_gen,
                                      n_mc,
                                      n_boot,
                                      progress_tf = TRUE) {

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
  if (n_neg > 0) {
    bootstrap_samples <- bootstrap_samples[bootstrap_samples[, 2] > 0, ]
    warning(paste0('Removing ', n_neg, ' bootstrap samples because sigma^* <= 0'), immediate. = TRUE)
    n_boot <- n_boot - n_neg
  }

  Lambda <- 1 + bootstrap_samples[, 3] *
    (thresh - bootstrap_samples[, 1])/bootstrap_samples[, 2]
  Lambda <- Lambda^(-1/bootstrap_samples[, 3])
  Lambda <- lt_gen*Lambda
  Lambda_orig <- 1 + k*(thresh - mu)/sigma
  Lambda_orig <- Lambda_orig^(-1/k)
  Lambda_orig <- lt_gen*Lambda_orig
  const <- Lambda/lt_gen

  bootstrap_samples <- bootstrap_samples[!is.na(Lambda), ]
  const <- const[!is.na(Lambda)]
  Lambda <- Lambda[!is.na(Lambda)]
  new_n_boot <- length(Lambda)
  if (new_n_boot < n_boot) {
    warning(paste0('Removing ', n_boot - n_new_boot,
                   ' bootstrap samples because Lambda^* is NA'), immediate. = TRUE)
    n_boot <- new_n_boot
  }

  prob_zero_obs <- dpois(0, Lambda)
  zero_prob_high <- prob_zero_obs > 0.9
  if (sum(zero_prob_high) > 0) {
    bootstrap_samples <- bootstrap_samples[!zero_prob_high, ]
    const <- const[!zero_prob_high]
    Lambda <- Lambda[!zero_prob_high]
    warning(paste0('Removing ', sum(zero_prob_high), ' bootstrap samples because the probability of zero threshold exceedances is > 90%'), immediate. = TRUE)
    n_boot <- n_boot - sum(zero_prob_high)
  }

  big_Lambda <- Lambda > 10*Lambda_orig
  if (sum(big_Lambda) > 0) {

    bootstrap_samples <- bootstrap_samples[!big_Lambda, ]
    Lambda <- Lambda[!big_Lambda]
    const <- const[!big_Lambda]
    warning(paste0('Removing ', sum(big_Lambda), ' bootstrap samples because Lambda* > 10Lambda_orig'), immediate. = TRUE)
    n_boot <- n_boot - sum(big_Lambda)
  }

  value <- list(par = c(mu, sigma, k),
                cov_mat = cov_mat,
                thresh = thresh,
                lt_gen = lt_gen,
                boot_samps = fullMaxDistUncertCpp(bootstrap_samples[, 1],
                                                  bootstrap_samples[, 2],
                                                  bootstrap_samples[, 3],
                                                  Lambda, const,
                                                  n_mc, n_boot, progress_tf))
  class(value) <- 'full_max_dist_uncert'
  value
}
