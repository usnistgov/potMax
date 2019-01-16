#'
#' @title Decluster a Time Series
#'
#' @description Decluster a time series by forming clusters and returning only
#'   cluster maximums
#'
#' @details Clusters are formed by sequential observations above the series mean
#'   value.  All observations below the mean value are discarded.  Cluster
#'   maximums are returned.
#'
#' @param complete_series (numeric vector) The time series.
#'
#' @param obs_times (numeric vector or NULL) If NULL, ignored; otherwise, the
#'   observed times of the cluster maximums are returned too.
#'
#' @return An S3 object of class \code{declustered_series} with components
#'
#'   \describe{
#'
#'   \item{\code{$declustered_series}}{The cluster maximums.}
#'
#'   \item{\code{$declustered_times}}{If \code{obs_times} is non NULL, the
#'   observed times of the cluster maximums; otherwise NULL.}
#'
#'   }
#'
#' @examples
#'
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
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
    if (length(complete_series) != length(obs_times)) {

      stop('must have length(complete_series) == length(obs_times)')
    }
  }

  series_mean <- mean(complete_series)
  y <- rep(NaN, length(complete_series))

  if (is.null(obs_times)) {

    declusterCpp(complete_series, y, series_mean)
    y <-  y[!is.nan(y)]
    value <- list(declustered_series = y,
                  declustered_times = NULL)
    class(value) <- 'declustered_series'
    value
  } else {

    if (length(complete_series) != length(obs_times)){
      stop('length(complete_series) should == length(obs_times)')
    }
    dt <- rep(NaN, length(obs_times))
    declusterWithTimeCpp(complete_series, obs_times,
                         y, dt, series_mean)
    y <- y[!is.nan(y)]
    dt <- dt[!is.nan(dt)]
    value <- list(declustered_series = y,
                  declustered_times = dt)
    class(value) <- 'declustered_series'
    value
  }
}
