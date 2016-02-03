# #' @title rGumbelPOT
# #'
# #' @description rGumbelPOT
# #'
# #' @details Generates a single POT series for the Gumbel (zero tail parameter)
# #'   intensity function
# #'
# #' @param mu Location parameter
# #'
# #' @param sigma Scale parameter
# #'
# #' @param thresh The threshold
# #'
# #' @param lt Lengt of the series in seconds
# #'
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

# #' @title rFullPOT
# #'
# #' @description rFullPOT
# #'
# #' @details Generates a single POT series for the full intensity function
# #'
# #' @param mu Location parameter
# #'
# #' @param sigma Scale parameter
# #'
# #' @param k Tail length parameter
# #'
# #' @param thresh The threshold
# #'
# #' @param lt Lengt of the series
# #'
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
