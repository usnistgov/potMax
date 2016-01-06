#' @title fullLogLike
#'
#' @description fullLogLike
#'
#' @details Defines the full log-likelihood for the Poisson process POT model
#'
#' @param theta The current value of the parameter triple
#'
#' @param y The vector of observations that exceed the threshold
#'
#' @param thresh The threshold
#'
#' @param lt - The length of time in seconds over which data were observed
#'
#' @param N - The number of observations (also the length of the y parameter)
#'
fullLogLike <- function (theta, y, thresh, lt, N) {

  ## unpack the parameters
  mu <- theta[1]
  sigma <- theta[2]
  k <- theta[3]

  ## missing or null parameter values return
  ## -Inf
  if (is.na(mu) || is.null(mu)) {

    return(-Inf)
  }
  if (is.na(sigma) || is.null(sigma)) {

    return(-Inf)
  }
  if (is.na(k) || is.null(k)) {

    return(-Inf)
  }

  # if the shape parameter is zero use limit of the likelihood
  if (round(k, 10) == 0) {

    # the scale parameter must be > zero
    if (sigma <= 0) {

      return(-Inf)
    } else {

      sum_y <- sum(y)

      term1 <- -N*log(sigma)

      term2 <- (-1/sigma)*sum_y

      term3 <- N*mu/sigma

      term4 <- -lt
      term4 <- term4*exp((-(thresh - mu))/sigma)

      value <- term1 + term2 + term3 + term4
    }
  } else {
    # if the shape parameter is not zero, use the regular likelihood

    # the scale parameter must be > zero
    if (sigma <= 0) {

      return(-Inf)
    } else {

      # some checks to make sure that the current parameters are consistent with
      # the threshold and data
      check <- 1 + k*((thresh - mu)/sigma)
      term2 <- 1 + k*((y - mu)/sigma)
      if (min(c(check, term2)) <= 0) {

        return(-Inf)
      } else {

        term2 <- sum(log(term2))
        term2 <- ((-1/k) - 1)*term2

        term1 <- (-1)*N*log(sigma)

        term3 <- 1 + k*((thresh - mu)/sigma)
        term3 <- term3^(-1/k)
        term3 <- lt*term3

        value <- term1 + term2 - term3
      }
    }
  }

  value
}

#'
#' @title fullMLE
#'
#' @description fullMLE
#'
#' @details  Maximizes the full Poisson process POT log-likelihood from Smith
#'   2004 (book chapter).
#'
#' @param y (numeric vector) The vector of observations that exceed the
#'   threshold
#'
#' @param thresh (numeric scalar) The threshold
#'
#' @param lt (numeric scalar) The length of time over which data were observed in
#'   seconds
#'
#' @param N - (numeric scalar) The number of observations that exceed the
#'   threshold (also the length of y)
#'
#' @param n_starts (numeric scalar) The number of random starts to use in the
#'   search for the maximum
#'
#' @param hessian (logical scalar) Compute the hessian or not
#'
#' @export
#'
fullMLE <- function (y, thresh, lt, N, n_starts, hessian) {

  mle <- NULL

  # Start at the maximum for the Gumbel model first

  start <- gumbelMLE(y = y, lt=lt, thresh=thresh, hessian = FALSE)$par
  start <- c(start, 0.0)
  tmp_mle <- try(optim(par=start, fn=fullLogLike,
                       control=list(fnscale=-1, maxit=10000),
                       hessian=FALSE, y=y, thresh=thresh,
                       lt=lt, N=N), FALSE)

  ## check for convergence of optim
  if (!inherits(tmp_mle, "try-error")) {

    if (tmp_mle$convergence == 0) {

      mle <- tmp_mle
    }
  }

  # next, we pass several random starts into optim, and we take the parameters
  # corresponding to the largest likelihood value.

  for (i in 1:n_starts) {

    ## find a random starting value that does not have
    ## negative infinity as its likelihood value
    repeat {

      tmp_start_mu <- rnorm(n=1, mean=start[1], sd=2)
      tmp_start_sigma <- abs(rnorm(n=1, mean=start[2], sd=2))
      tmp_start_k <- runif(n=1, min=-0.25, max=0.25)
      tmp_start <- c(tmp_start_mu,
                     tmp_start_sigma,
                     tmp_start_k)

      check <- fullLogLike(theta=tmp_start, y=y, thresh=thresh,
                           lt=lt, N=N)

      if (check > -Inf) {

        break()
      }
    }

    tmp_mle <- try(optim(par=tmp_start, fn=fullLogLike,
                         control=list(fnscale=-1, maxit=10000),
                         hessian=FALSE, y=y, thresh=thresh,
                         lt=lt, N=N), FALSE)

    ## check for convergence of optim and if
    ## the new parameter values lead to a larger
    ## likelihood than the old parameter values
    if (!inherits(tmp_mle, "try-error")) {

      if (tmp_mle$convergence == 0) {

        if (is.null(mle)) {

          mle <- tmp_mle
        } else {

          if (tmp_mle$value > mle$value) {

            mle <- tmp_mle
          }
        }
      }
    }
  }

  if (is.null(mle)) {

    warning('fullMLE unable to maximize the log-likelihood')
  }

  # if hessian is TRUE and optim converged we run optim one last time to get the
  # hessian matrix
  if (hessian & !is.null(mle)) {

    tmp_mle <- try(optim(par=mle$par, fn=fullLogLike,
                         control=list(fnscale=-1, maxit=10000),
                         hessian=TRUE, y=y, thresh=thresh,
                         lt=lt, N=N), FALSE)

    if (!inherits(tmp_mle, "try-error")) {

      if (tmp_mle$convergence == 0) {

        mle <- tmp_mle
      } else {

        mle <- NULL
        warning('fullMLE failed hessian computation')
      }
    } else {

      mle <- NULL
      warning('fullMLE failed hessian computation')
    }
  }

  mle
}
