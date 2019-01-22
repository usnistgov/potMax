#' @title Random Sample from the Gumbel like Poisson Process
#'
#' @description Generate a random sample from the 2D Poisson process with a
#'   Gumbel like intensity function
#'
#' @details The intensity function is described in \code{gumbelMLE}
#'
#' @param mu (numeric scalar) Location parameter
#'
#' @param sigma (numeric scalar) Scale parameter
#'
#' @param thresh (numeric scalar) The threshold
#'
#' @param lt (numeric scalar) Lengt of the series in units of time (seconds,
#'   minutes, hours, etc.)
#'
#' @return A list with components \code{$t} the times of the threshold crossings
#'   in the same units of \code{lt} and \code{$y} the value of the threshold
#'   crossings (NOT exceedances of/differences from the threshold)
#'
#' @examples
#'
#' \dontrun{
#'
#' g_pot_samp <- rGumbelPOT(mu = 50, sigma = 5, thresh = 55, lt = 100)
#'
#' }
#'
#' @export
#'
rGumbelPOT <- function (mu, sigma, thresh, lt) {

  Lambda <- lt*exp((-(thresh - mu))/sigma)

  n <- rpois(n = 1, lambda = Lambda)

  if (n == 0) {

    list(t = NULL, y = NULL)
  } else {

    t <- runif(n = n, min = 0, max = lt)
    z <- runif(n = n, min = 0, max = 1)
    const <- Lambda/lt
    y <- mu - sigma*(log(const) + log(1 - z))
    list(t = t, y = y)
  }
}

#' @title Random Sample from the Full Poisson Process
#'
#' @description Generate a random sample from the 2D Poisson process with the
#'   full intensity function
#'
#' @details The intensity function is described in \code{fullMLE}
#'
#' @param mu (numeric scalar) Location parameter
#'
#' @param sigma (numeric scalar) Scale parameter
#'
#' @param k (numeric scalar) Tail length parameter
#'
#' @param thresh (numeric scalar) The threshold
#'
#' @param lt (numeric scalar) Lengt of the series in units of time (seconds,
#'   minutes, hours, etc.)
#'
#' @return A list with components \code{$t} the times of the threshold crossings
#'   in the same units of \code{lt} and \code{$y} the value of the threshold
#'   crossings (NOT exceedances of/differences from the threshold)
#'
#' @examples
#'
#' \dontrun{
#'
#' f_pot_samp <- rFullPOT(mu = 50, sigma = 5, k = 0.2, thresh = 55, lt = 100)
#'
#' }
#'
#' @export
#'
rFullPOT <- function (mu, sigma, k, thresh, lt) {

  Lambda <- 1 + k*(thresh - mu)/sigma
  Lambda <- Lambda^(-1/k)
  Lambda <- lt*Lambda

  n <- rpois(n = 1, lambda = Lambda)

  if (n == 0) {

    list(t = NULL, y = NULL)
  } else {

    t <- runif(n = n, min = 0, max = lt)
    z <- runif(n = n, min = 0, max = 1)
    const <- Lambda/lt
    y <- (sigma/k)*((const*(1 - z))^(-k) - 1) + mu
    list(t = t, y = y)
  }
}
