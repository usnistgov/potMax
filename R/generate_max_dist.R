#' @title Distribution of the Maximum Using the Gumbel Model
#'
#' @description Empirically build the distribution of the maximum value over
#'   some user defined length of time assuming the underlying data generating
#'   mechanism is the 2D extremal Poisson process with the Gumbel like intensity
#'   function (i.e. the tail length is zero)
#'
#' @details The results of fitting a Gumbel like 2D extremal Poisson process are
#'   fed into this function.  Random processes are repeatedly generated from the
#'   fitted model, and the maximum of each random process is recorded.  The
#'   recoreded maximums represent an iid sample from the distribution of the
#'   maximum value for a process of the desired length.  Note that the desired
#'   length of the process can be different from the length of time over which
#'   the data used to fit the model were observed.
#'
#' @param x An S3 object of class \code{gumbel_pot_fit} or a numeric vector of
#'   lenght 2.  If the latter the first element of the vector should be the
#'   estimated location parameter \eqn{\mu}, and the second element should be the
#'   estimated scale parameter \eqn{\sigma}.
#'
#' @param lt_gen Length of each generated series.  The units (seconds, minutes,
#'   hours, etc.) should be consistent with the value of \code{lt} provided to
#'   \code{gumbelMLE}.
#'
#' @param n_mc The number of samples to draw from the distribution of the
#'   maximum
#'
#' @param progress_tf Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{gumbel_max_dist} with elements
#'
#' \describe{
#'   \item{\code{$par}}{The parameters used to generate the random processes}
#'
#'   \item{\code{$thres}}{The threshold used}
#'
#'   \item{\code{$lt_gen}}{The value of the \code{lt_gen} argument}
#'
#'   \item{\code{$max_dist}}{A numeric vector of length \code{n_mc} containing
#'   the samples from the distribution of the maximum}
#'
#' }
#' @examples
#'
#' \dontrun{
#'
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
#'
#' thresholded_obs <- gumbelEstThreshold(x = declustered_obs,
#'                                       lt = 100,
#'                                       n_min = 10,
#'                                       n_max = 100)
#'
#' gumbel_pot_fit <- gumbelMLE(x = thresholded_obs,
#'                             hessian_tf = TRUE)
#'
#' gumbel_max_dist <- gumbelMaxDist(x = gumbel_pot_fit, lt_gen = 200, n_mc = 1000)
#' }
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

#' @describeIn gumbelMaxDist
#'
#' @export
#'
gumbelMaxDist.gumbel_pot_fit <- function(x, lt_gen, n_mc,
                                         progress_tf = TRUE) {

  gumbelMaxDist.default(x = x$par,
                        thresh = x$thresh,
                        lt_gen = lt_gen,
                        n_mc = n_mc,
                        progress_tf = progress_tf)
}

#' @describeIn gumbelMaxDist
#'
#' @param thresh The threshold
#'
#' @export
#'
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

#' @title Uncertainty in the Distribution of the Maximum Using the Gumbel Model
#'
#' @description Evaluate uncertainty in the distribution of the maximum value
#'   over some user defined length of time assuming the underlying data
#'   generating mechanism is the 2D extremal Poisson process with the Gumbel
#'   like intensity function (i.e. the tail length is zero)
#'
#' @details The results of fitting a Gumbel like 2D extremal Poisson process are
#'   fed into this function.  The Hessian matrix is used to repeatedly perturb
#'   the MLE, and for each set of perturbed parameters the distribution of the
#'   maximum is empirically constructed as described in \code{gumbelMaxDist}.
#'   The bootstrap replicates of the distribution of the maximum may be used to
#'   quantify uncertainty and construct intervals.
#'
#' @param x An S3 object of class \code{gumbel_pot_fit} or a numeric vector of
#'   lenght 2.  If the latter the first element of the vector should be the
#'   estimated location parameter \eqn{\mu}, and the second element should be the
#'   estimated scale parameter \eqn{\sigma}.
#'
#' @param lt_gen Length of each generated series.  The units (seconds, minutes,
#'   hours, etc.) should be consistent with the value of \code{lt} provided to
#'   \code{gumbelMLE}.
#'
#' @param n_mc The number of samples to draw from the distribution of the
#'   maximum
#'
#' @param n_boot The number of bootstrap replicates of the distribution of the
#'   maximum to create.
#'
#' @param progress_tf Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{gumbel_max_dist_uncert} with elements
#'
#' \describe{
#'   \item{\code{$par}}{The parameters used to generate the random processes}
#'
#'   \item{\code{$cov_mat}}{The covariance matrix (negative inverse Hessian)
#'   used to perturb \code{$par}}
#'
#'   \item{\code{$thres}}{The threshold used}
#'
#'   \item{\code{$lt_gen}}{The value of the \code{lt_gen} argument}
#'
#'   \item{\code{$boot_samps}}{A matrix \code{n_boot} rows and \code{n_mc}
#'   columns containing the bootstrap replicates of the distribution of the
#'   maximum}
#'
#' }
#' @examples
#'
#' \dontrun{
#'
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
#'
#' thresholded_obs <- gumbelEstThreshold(x = declustered_obs,
#'                                       lt = 100,
#'                                       n_min = 10,
#'                                       n_max = 100)
#'
#' gumbel_pot_fit <- gumbelMLE(x = thresholded_obs,
#'                             hessian_tf = TRUE)
#'
#' gumbel_max_dist_uncert <- gumbelMaxDistUncert(x = gumbel_pot_fit, lt_gen = 200,
#'                                               n_mc = 1000, n_boot = 200)
#' }
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

#' @describeIn gumbelMaxDistUncert
#'
#' @export
#'
gumbelMaxDistUncert.gumbel_pot_fit <- function(x, lt_gen,
                                               n_mc,
                                               n_boot,
                                               progress_tf = TRUE) {

  gumbelMaxDistUncert.default(x = x$par,
                              cov_mat = -solve(x$lhessian),
                              thresh = x$thresh,
                              lt_gen = lt_gen,
                              n_mc = n_mc,
                              n_boot = n_boot,
                              progress_tf = progress_tf)
}

#' @describeIn gumbelMaxDistUncert
#'
#' @param cov_mat The covariance matrix to use to perturn the MLE (most usually
#'   the negative inverse of the Hessian matrix)
#'
#' @param thresh The threshold
#'
#' @export
#'
gumbelMaxDistUncert.default <- function(x, cov_mat,
                                        thresh, lt_gen,
                                        n_mc,
                                        n_boot,
                                        progress_tf = TRUE) {

  mu <- x[1]
  sigma <- x[2]
  lt_gen <- lt_gen[1]

  if (lt_gen <= 0) {
    stop('must have lt_gen > 0')
  }

  bootstrap_samples <- MASS::mvrnorm(n = n_boot, mu = c(mu, log(sigma)), Sigma = cov_mat)
  bootstrap_samples[, 2] <- exp(bootstrap_samples[, 2])

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

#' @title Distribution of the Maximum Using the Full Model
#'
#' @description Empirically build the distribution of the maximum value over
#'   some user defined length of time assuming the underlying data generating
#'   mechanism is the 2D extremal Poisson process with the full intensity
#'   function
#'
#' @details The results of fitting the 2D extremal Poisson process are
#'   fed into this function.  Random processes are repeatedly generated from the
#'   fitted model, and the maximum of each random process is recorded.  The
#'   recoreded maximums represent an iid sample from the distribution of the
#'   maximum value for a process of the desired length.  Note that the desired
#'   length of the process can be different from the length of time over which
#'   the data used to fit the model were observed.
#'
#' @param x An S3 object of class \code{full_pot_fit} or a numeric vector of
#'   lenght 3.  If the latter the first element of the vector should be the
#'   estimated location parameter \eqn{\mu}, the second element should be the
#'   estimated scale parameter \eqn{\sigma}, and the thrid element should be the
#'   estimated tail length parameter \eqn{k}.
#'
#' @param lt_gen Length of each generated series.  The units (seconds, minutes,
#'   hours, etc.) should be consistent with the value of \code{lt} provided to
#'   \code{gumbelMLE}.
#'
#' @param n_mc The number of samples to draw from the distribution of the
#'   maximum
#'
#' @param progress_tf Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{full_max_dist} with elements
#'
#' \describe{
#'   \item{\code{$par}}{The parameters used to generate the random processes}
#'
#'   \item{\code{$thres}}{The threshold used}
#'
#'   \item{\code{$lt_gen}}{The value of the \code{lt_gen} argument}
#'
#'   \item{\code{$max_dist}}{A numeric vector of length \code{n_mc} containing
#'   the samples from the distribution of the maximum}
#'
#' }
#' @examples
#'
#' \dontrun{
#'
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
#'
#' thresholded_obs <- fullEstThreshold(x = declustered_obs,
#'                                     lt = 100,
#'                                     n_min = 10,
#'                                     n_max = 100,
#'                                     n_starts = 10)
#'
#' full_pot_fit <- fullMLE(x = thresholded_obs,
#'                         hessian_tf = TRUE,
#'                         n_starts = 10)
#'
#' full_max_dist <- fullMaxDist(x = full_pot_fit, lt_gen = 200, n_mc = 1000)
#' }
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

#' @describeIn fullMaxDist
#'
#' @export
#'
fullMaxDist.full_pot_fit <- function(x, lt_gen, n_mc,
                                     progress_tf = TRUE) {

  fullMaxDist.default(x = x$par,
                      thresh = x$thresh,
                      lt_gen = lt_gen,
                      n_mc = n_mc,
                      progress_tf = progress_tf)
}

#' @describeIn fullMaxDist
#'
#' @param thresh The threshold
#'
#' @export
#'
fullMaxDist.default <- function(x, thresh,
                                lt_gen, n_mc,
                                progress_tf = TRUE) {

  mu <- x[1]
  sigma <- x[2]
  k <- x[3]
  lt_gen <- lt_gen[1] # only scalars are allowed

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

#' @title Uncertainty in the Distribution of the Maximum Using the Full Model
#'
#' @description Evaluate uncertainty in the distribution of the maximum value
#'   over some user defined length of time assuming the underlying data
#'   generating mechanism is the 2D extremal Poisson process with the Full
#'   intensity function
#'
#' @details The results of fitting the 2D extremal Poisson process are
#'   fed into this function.  The Hessian matrix is used to repeatedly perturb
#'   the MLE, and for each set of perturbed parameters the distribution of the
#'   maximum is empirically constructed as described in \code{fullMaxDist}.
#'   The bootstrap replicates of the distribution of the maximum may be used to
#'   quantify uncertainty and construct intervals.
#'
#' @param x An S3 object of class \code{full_pot_fit} or a numeric vector of
#'   lenght 3.  If the latter the first element of the vector should be the
#'   estimated location parameter \eqn{\mu}, the second element should be the
#'   estimated scale parameter \eqn{\sigma}, and the third element should be the
#'   estimated tail length parameter \eqn{k}.
#'
#' @param lt_gen Length of each generated series.  The units (seconds, minutes,
#'   hours, etc.) should be consistent with the value of \code{lt} provided to
#'   \code{gumbelMLE}.
#'
#' @param n_mc The number of samples to draw from the distribution of the
#'   maximum
#'
#' @param n_boot The number of bootstrap replicates of the distribution of the
#'   maximum to create.
#'
#' @param progress_tf Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{full_max_dist_uncert} with elements
#'
#' \describe{
#'   \item{\code{$par}}{The parameters used to generate the random processes}
#'
#'   \item{\code{$cov_mat}}{The covariance matrix (negative inverse Hessian)
#'   used to perturb \code{$par}}
#'
#'   \item{\code{$thres}}{The threshold used}
#'
#'   \item{\code{$lt_gen}}{The value of the \code{lt_gen} argument}
#'
#'   \item{\code{$boot_samps}}{A matrix \code{n_boot} rows and \code{n_mc}
#'   columns containing the bootstrap replicates of the distribution of the
#'   maximum}
#'
#' }
#' @examples
#'
#' \dontrun{
#'
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
#'
#' thresholded_obs <- fullEstThreshold(x = declustered_obs,
#'                                     lt = 100,
#'                                     n_min = 10,
#'                                     n_max = 100,
#'                                     n_starts = 10)
#'
#' full_pot_fit <- gumbelMLE(x = thresholded_obs,
#'                           hessian_tf = TRUE,
#'                           n_starts = 10)
#'
#' full_max_dist_uncert <- fullMaxDistUncert(x = full_pot_fit, lt_gen = 200,
#'                                           n_mc = 1000, n_boot = 200)
#' }
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

#' @describeIn fullMaxDistUncert
#'
#' @export
#'
fullMaxDistUncert.full_pot_fit <- function(x, lt_gen,
                                           n_mc,
                                           n_boot,
                                           progress_tf = TRUE) {

  fullMaxDistUncert.default(x = x$par,
                            cov_mat = -solve(x$lhessian),
                            thresh = x$thresh,
                            lt_gen = lt_gen,
                            n_mc = n_mc,
                            n_boot = n_boot,
                            progress_tf = progress_tf)
}

#' @describeIn fullMaxDistUncert
#'
#' @param cov_mat The covariance matrix to use to perturn the MLE (most usually
#'   the negative inverse of the Hessian matrix)
#'
#' @param thresh The threshold
#'
#' @export
#'
fullMaxDistUncert.default <- function(x,
                                      cov_mat, thresh, lt_gen,
                                      n_mc,
                                      n_boot,
                                      progress_tf = TRUE) {

  mu <- x[1]
  sigma <- x[2]
  k <- x[3]
  lt_gen <- lt_gen[1]

  if (lt_gen <= 0) {
    stop('must have lt_gen > 0')
  }

  bootstrap_samples <- MASS::mvrnorm(n = n_boot, mu = c(mu, log(sigma), k),
                                     Sigma = cov_mat)
  bootstrap_samples[, 2] <- exp(bootstrap_samples[, 2])

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
    warning(paste0('Removing ', n_boot - new_n_boot,
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
