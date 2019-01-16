#' @title Maximum Likelihood Estimation of the Gumble Model for Many Thresholds
#'
#' @description Fit the Gumbel like 2D extremal Poisson process for many
#'   thresholds
#'
#' @details \code{gumbelMLE} and \code{gumbelWPlot} are called for a sequence of
#'   thresholds.  Weights associated with each fit are also calculated.  Suppose
#'   that for threshold \eqn{u_i} the maximum vertical distance from a point on
#'   the W plot to the \eqn{45^\circ} line is \eqn{\delta_i} such that the
#'   \eqn{\delta_i} are scaled to the unit interval.  The weight
#'   associated with threshold \eqn{u_i} is then
#'
#'   \deqn{\frac{\exp\{-\tau\delta_i\}}{\sum\exp\{-\tau\delta_i\}}}
#'
#' @param x An S3 object of class \code{declustered_series} or a numeric vector.
#'   If the latter, the values to be thresholded and used in fitting.
#'
#' @param lt (numeric scalar) The length of the time series in units
#'   of time (seconds, minutes, hours, etc.).
#'
#' @param n_min (numeric scalar) The minimum number of thresholded observations
#'   to include
#'
#' @param n_max (numeric scalar) The maximum number of thresholded observations
#'   to include
#'
#' @param weight_scale (numeric scalar) The value of \eqn{\tau}
#'
#' @param progress_tf (logical scalar) Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{gumbel_multi_fit} with elements
#'
#' \describe{
#'   \item{\code{$all_fits}}{An object of type \code{gumbel_pot_fit} for each
#'   threshold}
#'
#'   \item{\code{$thresholds}}{The thresholds for the fits}
#'
#'   \item{\code{$weights}}{The weights associated with the fitted model for
#'   each threshold}
#'
#'   \item{\code{$lt}}{The value of the \code{lt} argument}
#'
#'   \item{\code{$n_min}}{The value of the \code{n_min} argument}
#'
#'   \item{\code{$n_max}}{The value of the \code{n_max} argument}
#'
#'   \item{\code{$weight_scale}}{The value of the \code{weight_scale} argument}
#' }
#'
#' @examples
#'
#' \dontrun{
#'
#' ddat <- decluster(-jp1tap813wind315$value)
#'
#' multi_est <- gumbelMultiFit(x = ddat, lt = 100, n_min = 10, n_max = 50, weight_scale = 5)
#'
#' }
#'
#' @export
#'
gumbelMultiFit <- function(x, lt, n_min, n_max, weight_scale,
                           progress_tf = TRUE) {
  UseMethod('gumbelMultiFit')
}

#' @describeIn gumbelMultiFit
#'
#' @export
#'
gumbelMultiFit.declustered_series <- function(x, lt,
                                              n_min,
                                              n_max,
                                              weight_scale,
                                              progress_tf = TRUE) {

  gumbelMultiFit.default(x = x$declustered_series,
                         lt = lt,
                         n_min = n_min,
                         n_max = n_max,
                         weight_scale = weight_scale,
                         progress_tf = progress_tf)
}

#' @describeIn gumbelMultiFit
#'
#' @export
#'
gumbelMultiFit.default <- function(x, lt, n_min, n_max, weight_scale,
                                   progress_tf = TRUE) {

  thresholds <- genThresholds(x, n_min, n_max)
  all_fits <- list()
  w_stats <- rep(NA, length(thresholds))

  if (progress_tf) {
    pb <- progress::progress_bar$new(total = length(thresholds), clear = FALSE,
                                     format = '|:bar| :percent ~ :eta',
                                     complete = '+', incomplete = ' ',
                                     current = ' ',
                                     width = floor(0.6*getOption('width')))
    pb$tick(0)
  }

  for (i in 1:length(thresholds)) {

    y <- x[x > thresholds[i]]

    all_fits[[i]] <- gumbelMLE.default(x = y,
                                       lt = lt,
                                       thresh = thresholds[i],
                                       hessian_tf = FALSE)
    w_stats[i] <- gumbelWPlot.default(x = all_fits[[i]]$par,
                                      y = y,
                                      thresh = thresholds[i],
                                      tf_plot = FALSE,
                                      BW = FALSE,
                                      details = FALSE)

    if (progress_tf) {
      pb$tick()
    }
  }

  min_w <- min(w_stats)
  max_w <- max(w_stats)
  tw_stats <- (weight_scale/(max_w - min_w))*(w_stats - min_w)
  tw_stats <- exp(-tw_stats)
  value <- list(all_fits = all_fits,
                w_stats = w_stats,
                thresholds = thresholds,
                weights = tw_stats/sum(tw_stats),
                lt = lt,
                n_min = n_min,
                n_max = n_max,
                weight_scale = weight_scale)
  class(value) <- 'gumbel_multi_fit'
  value
}

#' @title Maximum Likelihood Estimation of the Full Model for Many Thresholds
#'
#' @description Fit the full 2D extremal Poisson process for many thresholds
#'
#' @details \code{fullMLE} and \code{fullWPlot} are called for a sequence of
#'   thresholds.  Weights associated with each fit are also calculated.  Suppose
#'   that for threshold \eqn{u_i} the maximum vertical distance from a point on
#'   the W plot to the \eqn{45^\circ} line is \eqn{\delta_i} such that the
#'   \eqn{\delta_i} are scaled to the unit interval.  The weight
#'   associated with threshold \eqn{u_i} is then
#'
#'   \deqn{\frac{\exp\{-\tau\delta_i\}}{\sum\exp\{-\tau\delta_i\}}}
#'
#' @param x An S3 object of class \code{declustered_series} or a numeric vector.
#'   If the latter, the values to be thresholded and used in fitting.
#'
#' @param lt (numeric scalar) The length of the time series in units
#'   of time (seconds, minutes, hours, etc.).
#'
#' @param n_min (numeric scalar) The minimum number of thresholded observations
#'   to include
#'
#' @param n_max (numeric scalar) The maximum number of thresholded observations
#'   to include
#'
#' @param weight_scale (numeric scalar) The value of \eqn{\tau}
#'
#' @param n_starts (numeric scalar) The number of random starts to use in the
#'   search for the maximum
#'
#' @param progress_tf (logical scalar) Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{full_multi_fit} with elements
#'
#' \describe{
#'   \item{\code{$all_fits}}{An object of type \code{full_pot_fit} for each
#'   threshold}
#'
#'   \item{\code{$thresholds}}{The thresholds for the fits}
#'
#'   \item{\code{$weights}}{The weights associated with the fitted model for
#'   each threshold}
#'
#'   \item{\code{$lt}}{The value of the \code{lt} argument}
#'
#'   \item{\code{$n_min}}{The value of the \code{n_min} argument}
#'
#'   \item{\code{$n_max}}{The value of the \code{n_max} argument}
#'
#'   \item{\code{$weight_scale}}{The value of the \code{weight_scale} argument}
#' }
#'
#' @examples
#'
#' \dontrun{
#'
#' ddat <- decluster(-jp1tap813wind315$value)
#'
#' multi_est <- fullMultiFit(x = ddat, lt = 100, n_min = 10, n_max = 50, weight_scale = 5)
#'
#' }
#'
#' @export
#'
fullMultiFit <- function(x, lt, n_min, n_max, weight_scale, n_starts,
                         progress_tf) {
  UseMethod('fullMultiFit')
}

#' @describeIn fullMultiFit
#'
#' @export
#'
fullMultiFit.declustered_series <- function(x, lt,
                                            n_min,
                                            n_max,
                                            weight_scale,
                                            n_starts,
                                            progress_tf = TRUE) {

  fullMultiFit.default(x = x$declustered_series,
                       lt = lt,
                       n_min = n_min,
                       n_max = n_max,
                       weight_scale = weight_scale,
                       n_starts = n_starts,
                       progress_tf = progress_tf)
}

#' @describeIn fullMultiFit
#'
#' @export
#'
fullMultiFit.default <- function(x, lt, n_min, n_max,
                                 weight_scale, n_starts,
                                 progress_tf = TRUE) {

  thresholds <- genThresholds(x, n_min, n_max)
  all_fits <- list()
  w_stats <- rep(NA, length(thresholds))

  if (progress_tf) {
    pb <- progress::progress_bar$new(total = length(thresholds), clear = FALSE,
                                     format = '|:bar| :percent ~ :eta',
                                     complete = '+', incomplete = ' ',
                                     current = ' ',
                                     width = floor(0.6*getOption('width')))
    pb$tick(0)
  }

  for (i in 1:length(thresholds)) {

    y <- x[x > thresholds[i]]

    all_fits[[i]] <- fullMLE.default(x = y,
                                     lt = lt,
                                     thresh = thresholds[i],
                                     n_starts = n_starts,
                                     hessian_tf = FALSE)
    w_stats[i] <- fullWPlot.default(x = all_fits[[i]]$par,
                                    y = y,
                                    thresh = thresholds[i],
                                    tf_plot = FALSE,
                                    BW = FALSE,
                                    details = FALSE)

    if (progress_tf) {
      pb$tick()
    }
  }

  min_w <- min(w_stats)
  max_w <- max(w_stats)
  tw_stats <- weight_scale/(max_w - min_w)*(w_stats - min_w)
  tw_stats <- exp(-tw_stats)
  value <- list(all_fits = all_fits,
                w_stats = w_stats,
                thresholds = thresholds,
                weights = tw_stats/sum(tw_stats),
                lt = lt,
                n_min = n_min,
                n_max = n_max,
                weight_scale = weight_scale)
  class(value) <- 'full_multi_fit'
  value
}
