#' @title Return Values for the Gumbel Model
#'
#' @description The function calculates return values using the Gumbel like 2D
#'   extremal Poisson process
#'
#' @details Solves the equation
#'
#' \deqn{\int_{y_N}^\infty\int_0^1\lambda(t, y)dtdy = \frac{1}{N}}
#'
#' where \eqn{\lambda(t, y)} is given in the documentation for \code{gumbelMLE}
#'
#' @param x An S3 object of type \code{gumbel_pot_fit} or a numeric vector of
#'   length 2.  If the latter, the first element is interpreted as the location
#'   parameter \eqn{\mu}, and the second element is interpreted as the scale
#'   parameter \eqn{\sigma}.
#'
#' @param N (numeric scalar) The N in N-year return value.  This is a bit of a
#'   misnomer since the unit of time does not have to be years.  The function
#'   can calculate N-second, N-minute, N-hour, etc. return values as well.  In
#'   fact, the unit of time is the same unit of time passed in for the
#'   \code{lt} argument of \code{gumbelMLE}.  Naming the function for the unit
#'   of time year is simply due to my past experince with calculating return
#'   values on the time scale of years.
#'
#' @return An S3 object of class \code{gumbel_N_year_val} with elements
#'
#' \describe{
#'
#'   \item{\code{$par}}{numeric vector of length 2 containing the location and
#'   scale parameters, respectively of the 2D extremal Poisson process used to
#'   calculat the return value}
#'
#'   \item{\code{$thresh}}{The threshold}
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
#' thresholded_obs <- gumbelEstThreshold(x = declustered_obs,
#'                                       lt = 100,
#'                                       n_min = 10,
#'                                       n_max = 100)
#'
#' gumbel_pot_fit <- gumbelMLE(x = thresholded_obs,
#'                             hessian_tf = TRUE)
#'
#' 500_second_val <- gumbelNYear(x = gumbel_pot_fit, N = 500)
#' }
#'
#' @export
#'
gumbelNYear <- function(x, N, ...) {
  UseMethod('gumbelNYear')
}

#' @describeIn gumbelNYear
#'
#' @export
#'
gumbelNYear.gumbel_pot_fit <- function(x, N) {

  gumbelNYear.default(x = x$par,
                      thresh = x$thresh,
                      N = N)
}

#' @describeIn gumbelNYear
#'
#' @param thresh (numeric scalar) The threshold
#'
#' @export
#'
gumbelNYear.default <- function(x, thresh, N) {

  mu <- x[1]
  sigma <- x[2]

  N_year_val <- sigma*log(N) + mu

  value <- list(par = x,
                thresh = thresh,
                N = N,
                N_year_val = N_year_val)
  class(value) <- 'gumbel_N_year_val'
  value
}

#' @title Uncertainty in Return Values for the Gumbel Model
#'
#' @description The function produces a bootstrap sample of return values using
#'   the Gumbel like 2D extremal Poisson process
#'
#' @details Repeatedly solves the equation
#'
#' \deqn{\int_{y_N}^\infty\int_0^1\lambda(t, y)dtdy = \frac{1}{N}}
#'
#' where \eqn{\lambda(t, y)} is given in the documentation for \code{gumbelMLE}
#' for perturbed values of \eqn{\mu} and \eqn{\sigma}.  The perturbed values are
#' obtained by using the Hessian matrix and multivariate Guassian distribution
#' to perturb the MLE.
#'
#' @param x An S3 object of type \code{gumbel_pot_fit} or a numeric vector of
#'   length 2.  If the latter, the first element is interpreted as the location
#'   parameter \eqn{\mu}, and the second element is interpreted as the scale
#'   parameter \eqn{\sigma}.
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
#' @return An S3 object of class \code{gumbel_N_year_val_uncert} with elements
#'
#' \describe{
#'
#'   \item{\code{$par}}{numeric vector of length 2 containing the location and
#'   scale parameters, respectively of the 2D extremal Poisson process used to
#'   calculat the return value}
#'
#'   \item{\code{$thresh}}{The threshold}
#'
#'   \item{\code{$N}}{The value inputted for \code{N}}
#'
#'   \item{\code{$boot_samps}}{(numeric vector of length \code{n_boot}) The
#'   bootstrap sample of return values}
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
#' thresholded_obs <- gumbelEstThreshold(x = declustered_obs,
#'                                       lt = 100,
#'                                       n_min = 10,
#'                                       n_max = 100)
#'
#' gumbel_pot_fit <- gumbelMLE(x = thresholded_obs,
#'                             hessian_tf = TRUE)
#'
#' 500_second_val_uncert <- gumbelNYearUncert(x = gumbel_pot_fit, N = 500,
#'                                            n_boot = 1000)
#' }
#'
#' @export
#'
gumbelNYearUncert <- function(x, N, n_boot, ...) {
  UseMethod('gumbelNYearUncert')
}

#' @describeIn gumbelNYearUncert
#'
#' @export
#'
gumbelNYearUncert.gumbel_pot_fit <- function(x, N, n_boot) {

  gumbelNYearUncert.default(x = x$par,
                            cov_mat = -solve(x$lhessian),
                            thresh = x$thresh,
                            N = N,
                            n_boot = n_boot)
}

#' @describeIn gumbelNYearUncert
#'
#' @param cov_mat The covariance matrix used to perturb the estimated
#'   parameters.  This will most usually be the negative inverse of the Hessian
#'   matrix at the MLE
#'
#' @export
#'
gumbelNYearUncert.default <- function(x, cov_mat,
                                      thresh, N,
                                      n_boot) {

  mu <- x[1]
  sigma <- x[2]

  bootstrap_samples <- MASS::mvrnorm(n = n_boot, mu = c(mu, log(sigma)), Sigma = cov_mat)
  bootstrap_samples[, 2] <- exp(bootstrap_samples[, 2])

  N_year_vals <- bootstrap_samples[, 2]*log(N) + bootstrap_samples[, 1]

  value <- list(par = c(mu, sigma),
                cov_mat = cov_mat,
                thresh = thresh,
                N = N,
                boot_samps = N_year_vals)
  class(value) <- 'gumbel_N_year_val_uncert'
  value
}

# #' @title fullNYear
# #'
# #' @description fullNYear
# #'
# #' @details Calculates the N-year return value from a fitted POT model with
# #'   a non-zero tail length parameter
# #'
# #' @param x
# #'
# #' @param N The N in N-year return value
# #'
# #' @export
# #'
# #' @useDynLib potMax
# #'
# #' @importFrom Rcpp evalCpp
# #'
fullNYear <- function (x, N, ...) {
  UseMethod('fullNYear')
}

# #' @export
fullNYear.full_pot_fit <- function (x, N) {

  fullNYear.default(x = x$par,
                    N = N)
}

# #' @export
fullNYear.default <- function (x, N) {

  (((1/N)^(-x[3])) - 1)*(x[2]/x[3]) + x[1]
}
