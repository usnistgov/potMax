#'
#' @title gumbelMultiFit
#'
#' @description gumbelMultiFit
#'
#' @details Fit the Gumble POT model with multiple thresholds
#'
#' @param y_all The unthresholded series of observations
#'
#' @param lt The lenght of time in seconds over which the observations are
#'   recorded
#'
#' @param n_min The minimum number of thresholded observations to include
#'
#' @param n_max The maximum number of thresholded observations to include
#'
#' @param weight_scale  The weights are
#'   exp[-weight_scale/(max_w - min_w)*(w - min_w)]
#'   normalized to sum to unity
#'
#' @export
#'
gumbelMultiFit <- function(x, lt, n_min, n_max, weight_scale) {
  UseMethod('gumbelMultiFit')
}

#'
#' @export
#'
gumbelMultiFit.declustered_series <- function(x, lt,
                                              n_min,
                                              n_max,
                                              weight_scale) {

  gumbelMultiFit.default(x = x$declustered_series,
                         lt = lt,
                         n_min = n_min,
                         n_max = n_max,
                         weight_scale = weight_scale)
}

#' @export
gumbelMultiFit.default <- function(x, lt, n_min, n_max, weight_scale) {

  thresholds <- genThresholds(x, n_min, n_max)
  all_fits <- list()
  w_stats <- rep(NA, length(thresholds))

  pb <- progress::progress_bar$new(total = length(thresholds), clear = FALSE,
                                   complete = '*')
  pb$tick(0)

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

    pb$tick()
  }

  min_w <- min(w_stats)
  max_w <- max(w_stats)
  w_stats <- weight_scale/(max_w - min_w)*(w_stats - min_w)
  w_stats <- exp(-w_stats)
  value <- list(all_fits = all_fits,
                weights = w_stats/sum(w_stats),
                lt = lt,
                n_min = n_min,
                n_max = n_max,
                weight_scale = weight_scale)
  class(value) <- 'gumbel_multi_fit'
  value
}

#'
#' @title fullMultiFit
#'
#' @description fullMultiFit
#'
#' @details Fit the full POT model with multiple thresholds
#'
#' @param y_all The unthresholded series of observations
#'
#' @param lt The lenght of time in seconds over which the observations
#'   are recorded
#'
#' @param n_min The minimum number of thresholded observations to
#'   include
#'
#' @param n_max The maximum number of thresholded observations to
#'   include
#'
#' @param weight_scale  The weights are
#'   exp[-weight_scale/(max_w - min_w)*(w - min_w)]
#'   normalized to sum to unity
#'
#' @param n_starts (numeric scalar) The number of random starts to use
#'   in the search for the maximum
#'
#' @export
#'
fullMultiFit <- function(x, lt, n_min, n_max, weight_scale, n_starts) {
  UseMethod('fullMultiFit')
}

#'
#' @export
#'
fullMultiFit.declustered_series <- function(x, lt,
                                            n_min,
                                            n_max,
                                            weight_scale,
                                            n_starts) {

  fullMultiFit.default(x = x$declustered_series,
                       lt = lt,
                       n_min = n_min,
                       n_max = n_max,
                       weight_scale = weight_scale,
                       n_starts = n_starts)
}

#' @export
fullMultiFit.default <- function(x, lt, n_min, n_max,
                                 weight_scale, n_starts) {

  thresholds <- genThresholds(x, n_min, n_max)
  all_fits <- list()
  w_stats <- rep(NA, length(thresholds))

  pb <- progress::progress_bar$new(total = length(thresholds), clear = FALSE,
                                   complete = '*')
  pb$tick(0)

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

    pb$tick()
  }

  min_w <- min(w_stats)
  max_w <- max(w_stats)
  w_stats <- weight_scale/(max_w - min_w)*(w_stats - min_w)
  w_stats <- exp(-w_stats)
  value <- list(all_fits = all_fits,
                weights = w_stats/sum(w_stats))
  class(value) <- 'full_multi_fit'
  value
}