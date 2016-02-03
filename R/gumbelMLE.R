# #'
# #' @title sigmaScore
# #'
# #' @description sigmaScore
# #'
# #' @details This is one of the score equations for the the Poisson Process POT
# #'   Gumbel model (zero tail parameter), where the value of mu from solving the
# #'   other score equation, as a function of sigma is plugged in.  The root of
# #'   this function provides the MLE for sigma.
# #'
# #' @param sigma the value of the scale parameter at which the function is
# #'   evaluated
# #'
# #' @param N The length of the data vector
# #'
# #' @param lt The length of time in seconds over which observations were taken
# #'
# #' @param thresh The threshold over which all observations fall
# #'
# #' @param sum_y The sum of the observations
# #'
sigmaScore <- function (sigma, N, lt, thresh, sum_y) {

  if (sigma <= 0) {
    stop('sigma must be > 0')
  }

  if (N < 1) {
    stop('N must be >= 1')
  }

  if (lt <= 0) {
    stop('T must be > 0')
  }

  mu <- sigma*log((N/lt)) + thresh

  term1 <- -sigma*N

  term2 <- sum_y

  term3 <- -N*mu

  term4 <- -lt*(thresh - mu)
  term4 <- term4*exp((-(thresh - mu))/sigma)

  value <- term1 + term2 + term3 + term4
}

# #' @title sigmaMLE
# #'
# #' @description sigmaMLE
# #'
# #' @details Solves one of the two score equations for the limit of the Poisson
# #'   process likelihood where the value of mu from solving the other score
# #'   equation, as a function of sigma, is plugged in. The Poisson process
# #'   likelihood is actually the limit as the shape (tail) parameter goes to
# #'   zero.
# #'
# #' @param N The length of the data vector
# #'
# #' @param lt The length of time in seconds over which observations were taken
# #'
# #' @param thresh The threshold over which all observations fall
# #'
# #' @param sum_y The sum of the observations
# #'
sigmaMLE <- function (N, lt, thresh, sum_y) {

  if (N < 1) {
    stop('N must be >= 1')
  }

  if (lt <= 0) {
    stop('T must be > 0')
  }

  a <- 0.1

  b <- 1

  eq_a <- sigmaScore(sigma=a, lt=lt, N=N,
                     thresh=thresh, sum_y=sum_y)

  eq_b <- sigmaScore(sigma=b, lt=lt, N=N,
                     thresh=thresh, sum_y=sum_y)

  repeat {

    if(sign(eq_a) != sign(eq_b)) {

      break
    }

    a <- a/2

    eq_a <- sigmaScore(sigma=a, lt=lt, N=N,
                       thresh=thresh, sum_y=sum_y)

    if(sign(eq_a) != sign(eq_b)) {

      break
    }

    b <- b + 1

    eq_b <- sigmaScore(sigma=b, lt=lt, N=N,
                       thresh=thresh, sum_y=sum_y)

    if(sign(eq_a) != sign(eq_b)) {

      break
    }

    if (b > 100) {

      stop("A reasonable value for sigma does not exist")
    }
  }

  sigma <- uniroot(f=sigmaScore, interval=c(a, b),
                   lt=lt, N=N, thresh=thresh, sum_y=sum_y)
  sigma$root
}

#' @title Maximum Likelihood Estimation for the Gumble Model
#'
#' @description Solves the score equations for the Poisson process likelihood
#'   using the Gumbel intensity function
#'
#' @details  The likelihood is
#'
#'   \deqn{\left(\prod_{i = 1}^I \lambda(t_i,
#'   y_i)\right)\exp\left[-\int_\mathcal{D} \lambda(t, y)dtdy\right]}
#'
#'   where
#'
#'   \deqn{\lambda(t, y) = \frac{1}{\sigma}\exp\left[\frac{-(y -
#'   \mu)}{\sigma}\right]}
#'
#' @param x An S3 object of class \code{thresholded_series} or a numeric vector.
#'   If the latter, the values used in fitting.
#'
#' @param hessian_tf (logical scalar) Compute the Hessian matrix (TRUE) or not.
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
#' }
#'
#' @export
#'
gumbelMLE <- function (x, hessian_tf, ...) {
  UseMethod('gumbelMLE')
}

#' @describeIn gumbelMLE
#'
#' @export
#'
gumbelMLE.thresholded_series <- function (x, hessian_tf) {

  gumbelMLE.default(x = x$y,
                    lt = x$lt,
                    thresh = x$selected_threshold,
                    hessian_tf = hessian_tf)
}

#' @describeIn gumbelMLE
#'
#' @param lt (numeric scalar) The length of the time series in units
#'   of time (seconds, minutes, hours, etc.).
#'
#' @param thresh (numeric scalar) The threshold for the values
#'
#' @export
#'
gumbelMLE.default <- function (x, lt, thresh, hessian_tf) {

  N <- length(x)
  sigma <- sigmaMLE(lt=lt, N=N,
                    thresh=thresh, sum_y=sum(x))

  mu <- sigma*log((N/lt)) + thresh

  if (hessian_tf) {

    hess_est <- try(gumbelHessian(theta = c(mu, sigma),
                                  thresh = thresh, y = x,
                                  lt = lt, N = N))
    if (!inherits(hess_est, 'try-error')) {

      value <- list(par = c(mu, sigma), hessian = hess_est,
                    y = x, thresh = thresh)
      class(value) <- 'gumbel_pot_fit'
    } else {

      warning('Hessian computation failed')
      value <- list(par = c(mu, sigma), hessian = NULL,
                    y = x, thresh = thresh)
      class(value) <- 'gumbel_pot_fit'
    }
  } else {

    value <- list(par = c(mu, sigma), hessian = NULL,
                  y = x, thresh = thresh)
    class(value) <- 'gumbel_pot_fit'
  }

  value
}

# #'
# #' @title gumbelLogLike
# #'
# #' @description gumbelLogLike
# #'
# #' @details Defines the log-likelihood when the tail length parameter is exactly
# #'   zero
# #'
# #' @param theta The current value of the parameters (mu, sigma)
# #'
# #' @param y The vector of observations that exceed the threshold
# #'
# #' @param thresh The threshold
# #'
# #' @param lt The length of time over which data were observed in seconds
# #'
# #' @param N The number of observations that exceed the threshold (also the
# #'   length of y)
# #'
gumbelLogLike <- function (theta, y, thresh, lt, N) {

  full_theta <- c(theta, 0.0)

  fullLogLike(theta=full_theta,
              y=y,
              thresh=thresh,
              lt=lt,
              N=N)
}

# #'
# #' @title gumbelHessian
# #'
# #' @description gumbelHessian
# #'
# #' @details Calculates the Hessian matrix at the maximum value of the Poisson
# #'   process POT log-likelihood with the tail parameter fixed at zero
# #'
# #' @param theta The value of the parameter pair (mu, sigma) that maximizes the
# #'   Poisson process POT log-likelihood with the tail length parameter fixed at
# #'   zero
# #'
# #' @param y The vector of observations that exceed the threshold
# #'
# #' @param thresh The threshold
# #'
# #' @param lt The length of time over which data were observed in seconds
# #'
# #' @param N The number of observations that exceed the threshold (also the
# #'   length of y)
# #'
gumbelHessian <- function (theta, y, thresh, lt, N) {

  numDeriv::hessian(func=gumbelLogLike,
                    x=theta,
                    y=y,
                    thresh=thresh,
                    lt=lt,
                    N=N)
}
