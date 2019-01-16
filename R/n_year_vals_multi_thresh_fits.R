#' @title Return Values for the Gumbel Model Using Many Thresholds
#'
#' @description The function calculates return values using the Gumbel like 2D
#'   extremal Poisson process, but combining the results over many thresholds
#'
#' @details Solves the equation
#'
#' \deqn{\int_{y_N}^\infty\int_0^1\lambda(t, y)dtdy = \frac{1}{N}}
#'
#' where \eqn{\lambda(t, y)} is given in the documentation for \code{gumbelMLE}
#' for many thresholds and combines the results by weighted averaging.  The
#' weights are described in the documentation for \code{gumbelMultiFit}.
#'
#' @param x An S3 object of type \code{gumbel_multi_fit}
#'
#' @param N (numeric scalar) The N in N-year return value.  This is a bit of a
#'   misnomer since the unit of time does not have to be years.  The function
#'   can calculate N-second, N-minute, N-hour, etc. return values as well.  In
#'   fact, the unit of time is the same unit of time passed in for the \code{lt}
#'   argument of \code{gumbelMultiFit}.  Naming the function for the unit of
#'   time year is simply due to my past experince with calculating return values
#'   on the time scale of years.
#'
#' @return An S3 object of class \code{gumbel_N_year_val_multi_thresh} with
#'   elements
#'
#'
#' \describe{
#'
#'   \item{\code{$mu}}{A numeric vector containing the estimated location
#'   parameters, one for each threshold considered}
#'
#'   \item{\code{$sigma}}{A numeric vector containing the estimated scale
#'   parameters, one for each threshold considered}
#'
#'   \item{\code{$thresh}}{A numeric vector containing the thresholds being
#'   considered}
#'
#'   \item{\code{$N}}{The value inputted for \code{N}}
#'
#'   \item{\code{$N_year_val}}{(numeric scalar) The return value}
#' }
#'
#' @examples
#'
#' \dontrun{
#'
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
#'
#' gumbel_multi_fit <- gumbelMultiFit(x = declustered_obs,
#'                                    lt = 100,
#'                                    n_min = 10,
#'                                    n_max = 50,
#'                                    weight_scale = 5)
#'
#' 500_second_val <- gumbelNYear(x = gumbel_multi_fit, N = 500)
#' }
#'
#' @export
#'
gumbelNYear.gumbel_multi_fit <- function(x, N) {

  mu <- sapply(x$all_fits, function(x) x$par[1])
  sigma <- sapply(x$all_fits, function(x)x$par[2])
  thresh <- sapply(x$all_fits, function(x)x$thresh)
  N <- N[1] # only scalars are allowed

  if (sum(sigma <= 0) > 0) {
    stop('must have all sigma > 0')
  }
  if (N < 0) {
    stop('must have N >= 1')
  }

  N_year_val <- sigma*log(N) + mu
  N_year_val <- sum(N_year_val*x$weights)

  value <- list(mu = mu,
                sigma = sigma,
                thresh = thresh,
                N = N,
                N_year_val = N_year_val)
  class(value) <- 'gumbel_N_year_val_multi_thresh'
  value
}

#' @title Uncertainty in Return Values for the Gumbel Model Using Many
#'   Thresholds
#'
#' @description The function produces a bootstrap sample of return values using
#'   the Gumbel like 2D extremal Poisson process with many possible thresholds
#'
#' @details Repeatedly solves the equation
#'
#' \deqn{\int_{y_N}^\infty\int_0^1\lambda(t, y)dtdy = \frac{1}{N}}
#'
#' where \eqn{\lambda(t, y)} is given in the documentation for \code{gumbelMLE}
#' for many thresholds and \code{n_boot} bootstrap replicates of the
#' unthresholded data inputted to\code{gumbelMultiFit}.
#'
#' @param x An S3 object of type \code{gumbel_multi_fit}.
#'
#' @param declust_obs (numeric vector) The observed data used by
#'   \code{gumbelMultiFit}.  This will very often be the
#'   \code{$declustered_series} element of an S3 object of class
#'   \code{declustered_series}
#'
#' @param N (numeric scalar) The N in N-year return value.  This is a bit of a
#'   misnomer since the unit of time does not have to be years.  The function
#'   can calculate N-second, N-minute, N-hour, etc. return values as well.  In
#'   fact, the unit of time is the same unit of time passed in for the
#'   \code{lt} argument of \code{gumbelMLE}.  Naming the function for the unit
#'   of time year is simply due to my past experince with calculating return
#'   values on the time scale of years.
#'
#' @param n_boot (numeric scalar) The number of bootstrap samples
#'
#' @param progress_tf (logical scalar) Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{gumbel_N_year_val_uncert_multi_thresh}
#'   which is a list of length \code{n_boot} of S3 objects of class
#'   \code{gumbel_N_year_val_multi_thresh}.
#'
#' @examples
#'
#' \dontrun{
#'
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
#'
#' gumbel_multi_fit <- gumbelMultiFit(x = declustered_obs,
#'                                    lt = 100,
#'                                    n_min = 10,
#'                                    n_max = 50,
#'                                    weight_scale = 5)
#'
#' 500_second_val_uncert <- gumbelNYearUncert(x = gumbel_multi_fit,
#'                                            declust_obs = declustered_obs$declustered_series,
#'                                            N = 500,
#'                                            n_boot = 200)
#' }
#'
#' @export
#'
gumbelNYearUncert.gumbel_multi_fit <- function(x,
                                               declust_obs,
                                               N,
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

  boot_N_year <- list()

  for (i in 1:n_boot) {

    boot_multi_fit <- gumbelMultiFit.default(x = boot_samps[[i]],
                                             lt = x$lt,
                                             n_min = x$n_min,
                                             n_max = x$n_max,
                                             weight_scale = x$weight_scale,
                                             progress_tf = FALSE)
    boot_N_year[[i]] <- gumbelNYear(x = boot_multi_fit,
                                    N = N)
    if (progress_tf) {
      pb$tick()
    }
  }

  class(boot_N_year) <- 'gumbel_N_year_uncert_multi_thresh'
  boot_N_year
}
