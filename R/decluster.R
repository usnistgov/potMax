#'
#' @title decluster
#'
#' @description decluster
#'
#' @details Declusters a time series by taking maximums of clusters separated by
#'   mean value crossings.  Only clusters above the mean are counted.
#'
#' @param complete_series The time series
#'
#' @param obs_times
#'
#' @export
#'
#' @useDynLib potMax
#'
#' @importFrom Rcpp evalCpp
#'
decluster <- function (complete_series,
                       obs_times = NULL) {

  if (!is.null(obs_times)) {

    if (!is.numeric(complete_series) || !is.numeric((obs_times))) {

      stop('complete_series and obs_times must be numeric')
    }
  }

  series_mean <- mean(complete_series)
  y <- rep(NaN, length(complete_series))

  if (is.null(obs_times)) {

    declusterCpp(complete_series, y, series_mean)
    y <-  y[!is.nan(y)]
    list(declustered_series = y,
         declustered_times = NULL)
  } else {

    if (length(complete_series) != length(obs_times)){
      stop('length(complete_series) should == length(obs_times)')
    }
    dt <- rep(NaN, length(obs_times))
    declusterWithTimeCpp(complete_series, obs_times,
                         y, dt, series_mean)
    y <- y[!is.nan(y)]
    dt <- dt[!is.nan(dt)]
    list(declustered_series = y,
         declustered_times = dt)
  }
}
