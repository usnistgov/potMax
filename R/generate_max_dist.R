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
gumbelMaxDist <- function (mu, sigma, thresh, lt_gen, n_mc) {

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

  list(mu = mu,
       sigma = sigma,
       thresh = thresh,
       lt_gen = lt_gen,
       max_dist = gumbelMaxDistCpp(mu, sigma,
                                   Lambda, const, n_mc))
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
gumbelMaxDistUncert <- function (mu, sigma, cov_mat,
                                 thresh, lt_gen,
                                 n_mc, n_boot) {

  if (sigma <= 0) {
    stop('must have sigma > 0')
  }
  if (lt_gen <= 0) {
    stop('must have lt_gen > 0')
  }

  bootstrap_samples <- MASS::mvrnorm(n = n_boot, mu = c(mu, sigma), Sigma = cov_mat)
  if (sum(bootstrap_samples[, 2] <= 0) > 0){
    stop('negative bootstrap samples for sigma')
  }

  Lambda <- lt_gen*exp((-(thresh - bootstrap_samples[, 1]))/
                         bootstrap_samples[, 2])
  const <- Lambda/lt_gen

  prob_zero_obs <- dpois(0, Lambda)
  if (sum(prob_zero_obs > 0.9) > 0) {

    warning('The probability of zero observations for some of the bootstrap samples of the parameters is greater than 0.9')
  }

  list(mu = mu,
       sigma = sigma,
       cov_mat = cov_mat,
       thresh = thresh,
       lt_gen = lt_gen,
       boot_samps = gumbelMaxDistUncertCpp(bootstrap_samples[, 1],
                                           bootstrap_samples[, 2],
                                           Lambda, const,
                                           n_mc, n_boot))
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
fullMaxDist <- function (mu, sigma, k, thresh,
                         lt_gen, n_mc) {

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

  list(mu = mu,
       sigma = sigma,
       k = k,
       thresh = thresh,
       lt_gen = lt_gen,
       max_dist = fullMaxDistCpp(mu, sigma, k,
                                 Lambda, const, n_mc))
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
fullMaxDistUncert <- function (mu, sigma, k,
                               cov_mat, thresh, lt_gen,
                               n_mc, n_boot) {

  if (sigma <= 0) {
    stop('must have sigma > 0')
  }
  if (lt_gen <= 0) {
    stop('must have lt_gen > 0')
  }

  bootstrap_samples <- MASS::mvrnorm(n = n_boot, mu = c(mu, sigma, k),
                               Sigma = cov_mat)
  if (sum(bootstrap_samples[, 2] <= 0) > 0){
    stop('negative bootstrap samples for sigma')
  }

  Lambda <- 1 + bootstrap_samples[, 3]*
    (thresh - bootstrap_samples[, 1])/bootstrap_samples[, 2]
  Lambda <- Lambda^(-1/bootstrap_samples[, 3])
  Lambda <- lt_gen*Lambda
  const <- Lambda/lt_gen

  prob_zero_obs <- dpois(0, Lambda)
  if (sum(prob_zero_obs > 0.9) > 0) {

    warning('The probability of zero observations for some of the bootstrap samples of the parameters is greater than 0.9')
  }

  list(mu = mu,
       sigma = sigma,
       k = k,
       cov_mat = cov_mat,
       thresh = thresh,
       lt_gen = lt_gen,
       boot_samps = fullMaxDistUncertCpp(bootstrap_samples[, 1],
                                         bootstrap_samples[, 2],
                                         bootstrap_samples[, 3],
                                         Lambda, const,
                                         n_mc, n_boot))
}
