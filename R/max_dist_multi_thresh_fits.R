#' @title Distribution of the Maximum Using the Gumbel Model with Many Thresholds
#'
#' @description Empirically build the distribution of the maximum value over
#'   some user defined length of time assuming the underlying data generating
#'   mechanism is a 2D extremal Poisson process with the Gumbel like intensity
#'   function (i.e. the tail length is zero), but without specifying a single
#'   threshold.
#'
#' @details The results of fitting a Gumbel like 2D extremal Poisson process for
#'   many thresholds are fed into this function.  Random processes are
#'   repeatedly generated from each fitted model, according to the weights
#'   described in \code{gumbelMultiFit}, and the maximum of each random process
#'   is recorded.  The recoreded maximums represent an iid sample from a mixture
#'   of the distributions of maximum values.  Each potential threshold gives
#'   rise to a different distribution of the maximum value. Note that the
#'   desired length of the processes generated can be different from the length
#'   of time over which the data used to fit the models were observed.
#'
#' @param x An S3 object of class \code{gumbel_multi_fit}.
#'
#' @param lt_gen (numeric scalar) Length of each generated series.  The units
#'   (seconds, minutes, hours, etc.) should be consistent with the value of
#'   \code{lt} provided to \code{gumbelMultiFit}.
#'
#' @param n_mc (numeric scalar) The number of samples to draw from the mixture
#'   distribution of the maximums
#'
#' @param progress_tf (logical scalar) Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{gumbel_max_dist_multi_thresh} with elements
#'
#' \describe{
#'   \item{\code{$mu}}{The estimated location parameter for each threshold}
#'
#'   \item{\code{$sigma}}{The estimated scale parameter for each threshold}
#'
#'   \item{\code{$thres}}{The thresholds}
#'
#'   \item{\code{$lt_gen}}{The value of the \code{lt_gen} argument}
#'
#'   \item{\code{$n_each}}{The number of samples drawn from each mixture
#'   component.  The sum of \code{$n_each} should be \code{n_mc}}
#'
#'   \item{\code{$max_dist}}{A numeric vector of length \code{n_mc} containing
#'   the samples from the mixture distribution of the maximums}
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
#' gumbel_multi_fit <- gumbelMultiFit(x = declustered_obs, lt = 100,
#'                                    n_min = 10, n_max = 50,
#'                                    weight_scale = 5)
#'
#' gumbel_max_dist <- gumbelMaxDist(x = gumbel_multi_fit, lt_gen = 200, n_mc = 1000)
#' }
#'
#' @export
#'
#' @useDynLib potMax
#'
#' @importFrom Rcpp evalCpp
#'
gumbelMaxDist.gumbel_multi_fit <- function(x,
                                           lt_gen, n_mc,
                                           progress_tf = TRUE) {

  mu <- sapply(x$all_fits, function(x) x$par[1])
  sigma <- sapply(x$all_fits, function(x)x$par[2])
  thresh <- sapply(x$all_fits, function(x)x$thresh)
  lt_gen <- lt_gen[1] # only scalars are allowed

  if (sum(sigma <= 0) > 0) {
    stop('must have all sigma > 0')
  }
  if (lt_gen <= 0) {
    stop('must have lt_gen > 0')
  }

  Lambda <- lt_gen*exp((-(thresh - mu))/sigma)
  const <- Lambda/lt_gen

  prob_zero_obs <- dpois(0, Lambda)
  if (sum(prob_zero_obs > 0.9) > 0) {

    warning(paste0('The probability of zero observations is more than 0.9 for ', sum(prob_zero_obs > 0.9), ' thresholds\n'), immediate. = TRUE)
  }

  n_each <- round(n_mc*x$weights)
  long_mu <- rep(mu, n_each)
  long_sigma <- rep(sigma, n_each)
  long_thresh <- rep(thresh, n_each)
  long_Lambda <- rep(Lambda, n_each)
  long_const <- rep(const, n_each)

  value <- list(mu = mu,
                sigma = sigma,
                thresh = thresh,
                lt_gen = lt_gen,
                n_each = n_each,
                max_dist = gumbelMaxDistMultiCpp(long_mu, long_sigma,
                                                 long_Lambda,
                                                 long_const,
                                                 sum(n_each),
                                                 progress_tf))
  class(value) <- 'gumbel_max_dist_multi_thresh'
  value
}

#' @title Uncertainty in the Distribution of the Maximum Using the Gumbel Model
#'   With Many Thresholds
#'
#' @description Evaluate uncertainty in the mixture of distributions of maximums
#'   assuming the underlying data generating mechanism is the 2D extremal
#'   Poisson process with the Gumbel like intensity function (i.e. the tail
#'   length is zero), but without specifying a single threshold.  The mixture is
#'   across multiple thresholds.
#'
#' @details The results of fitting a many Gumbel like 2D extremal Poisson
#'   process are fed into this function.  The declustered data are repeatedly
#'   sampled with replacement, and for each resampled data set the distribution
#'   of the maximum is empirically constructed as described in
#'   \code{gumbelMaxDist.gumbel_multi_fit}. The bootstrap replicates of the
#'   mixture of distributions of the maximum may be used to quantify uncertainty
#'   and construct intervals.
#'
#' @param x An S3 object of class \code{gumbel_multi_fit}.
#'
#' @param declust_obs (numeric vector) The observed data used by
#'   \code{gumbelMultiFit}.  This will very often be the
#'   \code{$declustered_series} element of an S3 object of class
#'   \code{declustered_series}
#'
#' @param lt_gen (numeric scalar) Length of each generated series.  The units
#'   (seconds, minutes, hours, etc.) should be consistent with the value of
#'   \code{lt} provided to \code{gumbelMLE}.
#'
#' @param n_mc (numeric scalar) The number of samples to draw from the
#'   distribution of the maximum
#'
#' @param n_boot (numeric scalar) The number of bootstrap replicates of the
#'   distribution of the maximum to create.
#'
#' @param progress_tf Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{gumbel_max_dist_uncert_multi_thresh},
#'   which is a list of length \code{n_boot} of S3 objects of class
#'   \code{gumbel_max_dist_multi_thresh}.
#'
#' @examples
#'
#' \dontrun{
#'
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
#'
#' gumbel_multi_fit <- gumbelMultiFit(x = declustered_obs, lt = 100,
#'                                    n_min = 10, n_max = 50,
#'                                    weight_scale = 5)
#'
#' gumbel_multi_fit_uncert <- gumbelMaxDistUncert(x = gumbel_multi_fit,
#'                                                declust_obs = declustered_obs$declustered_series,
#'                                                lt_gen = 200,
#'                                                n_mc = 1000,
#'                                                n_boot = 200)
#' }
#'
#' @export
#'
gumbelMaxDistUncert.gumbel_multi_fit <- function(x,
                                                 declust_obs,
                                                 lt_gen,
                                                 n_mc,
                                                 n_boot,
                                                 progress_tf = TRUE) {

  if (progress_tf) {
    pb <- progress::progress_bar$new(total = n_boot, clear = FALSE,
                                     format = '|:bar| :percent ~ :eta',
                                     complete = '+', incomplete = ' ',
                                     current = ' ',
                                     width = floor(0.6*getOption('width')))
    pb$tick(0)
  }

  boot_samps <- lapply(1:n_boot, function(i, x)sample(x,length(x),TRUE),
                       x = declust_obs)

  boot_max_dist <- list()

  for (i in 1:n_boot) {

    boot_multi_fit <- gumbelMultiFit.default(x = boot_samps[[i]],
                                             lt = x$lt,
                                             n_min = x$n_min,
                                             n_max = x$n_max,
                                             weight_scale = x$weight_scale,
                                             progress_tf = FALSE)
    boot_max_dist[[i]] <- gumbelMaxDist(x = boot_multi_fit,
                                        lt_gen = lt_gen, n_mc = n_mc,
                                        progress_tf = FALSE)
    if (progress_tf) {
      pb$tick()
    }
  }

  class(boot_max_dist) <- 'gumbel_max_dist_uncert_multi_thresh'
  boot_max_dist
}

#' @title Distribution of the Maximum Using the Full Model with Many Thresholds
#'
#' @description Empirically build the distribution of the maximum value over
#'   some user defined length of time assuming the underlying data generating
#'   mechanism is a 2D extremal Poisson process, but without specifying a single
#'   threshold.
#'
#' @details The results of fitting a 2D extremal Poisson process for many
#'   thresholds are fed into this function.  Random processes are repeatedly
#'   generated from each fitted model, according to the weights described in
#'   \code{fullMultiFit}, and the maximum of each random process is recorded.
#'   The recoreded maximums represent an iid sample from a mixture of the
#'   distributions of maximum values.  Each threshold gives rise to a different
#'   distribution of the maximum value. Note that the desired length of the
#'   processes generated can be different from the length of time over which the
#'   data used to fit the models were observed.
#'
#' @param x An S3 object of class \code{full_multi_fit}.
#'
#' @param lt_gen (numeric scalar) Length of each generated series.  The units
#'   (seconds, minutes, hours, etc.) should be consistent with the value of
#'   \code{lt} provided to \code{fullMultiFit}.
#'
#' @param n_mc (numeric scalar) The number of samples to draw from the mixture
#'   distribution of the maximums
#'
#' @param progress_tf (logical scalar) Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{full_max_dist_multi_thresh} with elements
#'
#' \describe{
#'   \item{\code{$mu}}{The estimated location parameter for each threshold}
#'
#'   \item{\code{$sigma}}{The estimated scale parameter for each threshold}
#'
#'   \item{\code{$k}}{The estimated tail length parameter for each threshold}
#'
#'   \item{\code{$thres}}{The thresholds}
#'
#'   \item{\code{$lt_gen}}{The value of the \code{lt_gen} argument}
#'
#'   \item{\code{$n_each}}{The number of samples drawn from each mixture
#'   component.  The sum of \code{$n_each} should be \code{n_mc}}
#'
#'   \item{\code{$max_dist}}{A numeric vector of length \code{n_mc} containing
#'   the samples from the mixture distribution of the maximums}
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
#' full_multi_fit <- fullMultiFit(x = declustered_obs, lt = 100,
#'                                n_min = 10, n_max = 50,
#'                                weight_scale = 5, n_starts = 20)
#'
#' full_max_dist <- fullMaxDist(x = full_multi_fit, lt_gen = 200, n_mc = 1000)
#' }
#'
#' @export
#'
#' @useDynLib potMax
#'
#' @importFrom Rcpp evalCpp
#'
fullMaxDist.full_multi_fit <- function(x,
                                       lt_gen, n_mc,
                                       progress_tf = TRUE) {

  mu <- sapply(x$all_fits, function(x) x$par[1])
  sigma <- sapply(x$all_fits, function(x)x$par[2])
  k <- sapply(x$all_fits, function(x)x$par[3])
  thresh <- sapply(x$all_fits, function(x)x$thresh)
  lt_gen <- lt_gen[1] # only scalars are allowed

  if (sum(sigma <= 0) > 0) {
    stop('must have all sigma > 0')
  }
  if (lt_gen <= 0) {
    stop('must have lt_gen > 0')
  }

  Lambda <- 1 + k*(thresh - mu)/sigma
  Lambda <- Lambda^(-1/k)
  Lambda <- lt_gen*Lambda
  const <- Lambda/lt_gen

  prob_zero_obs <- dpois(0, Lambda)
  if (sum(prob_zero_obs > 0.9) > 0) {

    warning(paste0('The probability of zero observations is more than 0.9 for ', sum(prob_zero_obs > 0.9), ' thresholds\n'), immediate. = TRUE)
  }

  n_each <- round(n_mc*x$weights)
  long_mu <- rep(mu, n_each)
  long_sigma <- rep(sigma, n_each)
  long_k <- rep(k, n_each)
  long_thresh <- rep(thresh, n_each)
  long_Lambda <- rep(Lambda, n_each)
  long_const <- rep(const, n_each)

  value <- list(mu = mu,
                sigma = sigma,
                k = k,
                thresh = thresh,
                lt_gen = lt_gen,
                n_each = n_each,
                max_dist = fullMaxDistMultiCpp(long_mu, long_sigma, long_k,
                                               long_Lambda,
                                               long_const,
                                               sum(n_each),
                                               progress_tf))
  class(value) <- 'full_max_dist_multi_thresh'
  value
}

#' @title Uncertainty in the Distribution of the Maximum Using the Full Model
#'   With Many Thresholds
#'
#' @description Evaluate uncertainty in the mixture of distributions of maximums
#'   assuming the underlying data generating mechanism is the 2D extremal
#'   Poisson process, but without specifying a single threshold.  The mixture is
#'   across multiple thresholds.
#'
#' @details The results of fitting many 2D extremal Poisson process are fed into
#'   this function.  The declustered data are repeatedly sampled with
#'   replacement, and for each resampled data set the distribution of the
#'   maximum is empirically constructed as described in
#'   \code{fullMaxDist.full_multi_fit}. The bootstrap replicates of the mixture
#'   of distributions of the maximum may be used to quantify uncertainty and
#'   construct intervals.
#'
#' @param x An S3 object of class \code{full_multi_fit}.
#'
#' @param declust_obs (numeric vector) The observed data used by
#'   \code{fullMultiFit}.  This will very often be the
#'   \code{$declustered_series} element of an S3 object of class
#'   \code{declustered_series}
#'
#' @param lt_gen (numeric scalar) Length of each generated series.  The units
#'   (seconds, minutes, hours, etc.) should be consistent with the value of
#'   \code{lt} provided to
#'   \code{gumbelMLE}.
#'
#' @param n_mc (numeric scalar) The number of samples to draw from the
#'   distribution of the maximum
#'
#' @param n_boot (numeric scalar) The number of bootstrap replicates of the
#'   distribution of the maximum to create.
#'
#' @param n_starts (numeric scalar) The number of random starts to use in the
#'   search for the maximum
#'
#' @param progress_tf (logical scalar) Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{full_max_dist_uncert_multi_thresh},
#'   which is a list of length \code{n_boot} of S3 objects of class
#'   \code{full_max_dist_multi_thresh}.
#'
#' @examples
#'
#' \dontrun{
#'
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
#'
#' full_multi_fit <- fullMultiFit(x = declustered_obs, lt = 100,
#'                                n_min = 10, n_max = 50,
#'                                weight_scale = 5, n_starts = 20)
#'
#' full_multi_fit_uncert <- fullMaxDistUncert(x = full_multi_fit,
#'                                            declust_obs = declustered_obs$declustered_series,
#'                                            lt_gen = 200,
#'                                            n_mc = 1000,
#'                                            n_boot = 200)
#' }
#'
#' @export
#'
fullMaxDistUncert.full_multi_fit <- function(x,
                                             declust_obs,
                                             lt_gen,
                                             n_mc,
                                             n_boot,
                                             n_starts,
                                             progress_tf = TRUE) {

  boot_samps <- lapply(1:n_boot, function(i, x)sample(x,length(x),TRUE),
                       x = declust_obs)
  if (progress_tf) {
    pb <- progress::progress_bar$new(total = n_boot, clear = FALSE,
                                     format = '|:bar| :percent ~ :eta',
                                     complete = '+', incomplete = ' ',
                                     current = ' ',
                                     width = floor(0.6*getOption('width')))
    pb$tick(0)
  }

  boot_samps <- lapply(1:n_boot, function(i, x)sample(x,length(x),TRUE),
                       x = declust_obs)

  boot_max_dist <- list()

  for (i in 1:n_boot) {

    boot_multi_fit <- fullMultiFit.default(x = boot_samps[[i]],
                                           lt = x$lt,
                                           n_min = x$n_min,
                                           n_max = x$n_max,
                                           weight_scale = x$weight_scale,
                                           n_starts = n_starts,
                                           progress_tf = FALSE)
    boot_max_dist[[i]] <- fullMaxDist(x = boot_multi_fit,
                                      lt_gen = lt_gen, n_mc = n_mc,
                                      progress_tf = FALSE)
    if (progress_tf) {
      pb$tick()
    }
  }

  class(boot_max_dist) <- 'full_max_dist_uncert_multi_thresh'
  boot_max_dist
}