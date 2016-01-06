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
#' @export
#'
#' @useDynLib potMax
#'
#' @importFrom Rcpp evalCpp
#'
decluster <- function (complete_series) {

  series_mean <- mean(complete_series)
  y <- rep(NaN, length(complete_series))

  declusterCpp(complete_series, y, series_mean)

  y[!is.nan(y)]
}
