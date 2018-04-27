#' @title fullWPlot
#'
#' @description fullWPlot
#'
#' @details this is the W plot from the 2004 Richard Smith book chapter
#'
#' @param mu
#'
#' @param sigma
#'
#' @param k
#'
#' @param y
#'
#' @param thresh
#'
#' @param tf_plot
#'
#' @param BW (logical) only use black and white for the plot.  Ignored if
#'   tf_plot == FALSE
#'
#' @param details
#'
#' @param ...
#'
#' @export
#'
fullWPlot <- function (x, tf_plot, BW, details, ...) {
  UseMethod('fullWPlot')
}

#' @export
fullWPlot.full_pot_fit <- function (x,
                                    tf_plot,
                                    BW,
                                    details,
                                    ...) {
  invisible(fullWPlot.default(x = x$par,
                              y = x$y,
                              thresh = x$thresh,
                              tf_plot = tf_plot,
                              BW = BW,
                              details = details,
                              ...))
}

#' @export
fullWPlot.default <- function(x,
                              y,
                              thresh,
                              tf_plot,
                              BW,
                              details,
                              ...) {

  mu <- x[1]
  sigma <- x[2]
  k <- x[3]

  # if the absolute value of a number is less than this consider it zero
  eps <- 1e-8

  # vectorize everything so that loops are unneccesary
  n <- length(y)
  if (n <= 0) {
    stop('No observations!')
  }
  vec_mu <- rep(mu, n)
  vec_sigma <- rep(sigma, n)
  vec_k <- rep(k, n)

  # the formula requires the excesses instead of the actual values
  excesses <- y - thresh

  ## calculate the W statistics
  if (abs(k) > eps) {

    W1 <- (vec_k*excesses)/(vec_sigma + vec_k*(thresh - vec_mu))
    W2 <- 1 + W1
    ## in case any part of W is less than 0 at this point,
    ## we raise it to 0
    W3 <- W2
    W3[W3 < 0] <- 0
    W4 <- log(W3)
    W <- (1/vec_k)*W4
  } else {

    W1 <- rep(NA, n)
    W2 <- rep(NA, n)
    W3 <- rep(NA, n)
    W4 <- rep(NA, n)
    W <- excesses/vec_sigma
  }

  ## get the values ready to return if details
  ## are asked for
  if (details) {


    y <- y[order(W)]
    excesses <- excesses[order(W)]
    W1 <- W1[order(W)]
    W2 <- W2[order(W)]
    W3 <- W3[order(W)]
    W4 <- W4[order(W)]
  }
  W <- sort(x = W)

  ## calculate the appropriate exp(1)
  ## quantiles
  quantiles <- ((1:n) - 0.375)/(n + 0.25)
  exp1_quantiles <- qexp(p = quantiles)

  if (tf_plot) {

    plot(x = exp1_quantiles, y = W,
         xlab = "exp(1) quantiles",
         ylab = "Ordered W-Statistics",
         main = "W-Statistic Plot",
         bty = "l",
         ...)

    if (BW) {

      abline(a = 0, b = 1, col = "black")
    } else {

      abline(a = 0, b = 1, col = "red")
    }
  }

  if (details) {

    value <- abs((W - exp1_quantiles))
    data.frame(y, excesses, W1, W2, W3, W4, W, value)
  } else {

    invisible(max(abs((W - exp1_quantiles))))
  }
}

#' @title gumbelWPlot
#'
#' @description gumbelWPlot
#'
#' @details this is the W plot from the 2004 Richard Smith book chapter
#'
#' @param mu
#'
#' @param sigma
#'
#' @param y
#'
#' @param thresh
#'
#' @param tf_plot
#'
#' @param BW (logical) only use black and white for the plot.  Ignored if
#'   tf_plot == FALSE
#'
#' @param details
#'
#' @param ...
#'
#' @export
#'
gumbelWPlot <- function (x, tf_plot, BW, details, ...) {
  UseMethod('gumbelWPlot')
}

#' @export
gumbelWPlot.gumbel_pot_fit <- function (x,
                                        tf_plot,
                                        BW,
                                        details,
                                        ...) {

  invisible(gumbelWPlot.default(x = x$par,
                                y = x$y,
                                thresh = x$thresh,
                                tf_plot = tf_plot,
                                BW = BW,
                                details = details,
                                ...))
}

#' @export
gumbelWPlot.default <- function (x,
                                 y,
                                 thresh,
                                 tf_plot,
                                 BW,
                                 details,
                                 ...) {

  invisible(fullWPlot.default(x = c(x, 0),
                              y = y,
                              thresh = thresh,
                              tf_plot = tf_plot,
                              BW = BW,
                              details = details,
                              ...))
}

# #' @title genThresholds
# #'
# #' @description genThresholds
# #'
# #' @details Generate the thresholds from which the estimated threshold will be
# #'   chosen
# #'
# #' @param y_all
# #'
# #' @param n_min
# #'
# #' @param n_max
# #'
genThresholds <- function(y_all, n_min, n_max) {

  y <- sort(y_all, decreasing = TRUE)
  thresholds <- NULL

  n_max <- min(n_max, (length(y) - 1))

  y1 <- y[n_min:n_max]
  y2 <- y[(n_min + 1):(n_max + 1)]

  thresholds <- unique((y1 + y2)/2)
}

#'
#' @title gumbelEstThreshold
#'
#' @description gumbelEstThreshold
#'
#' @details Select a threshold when using the Gumble POT modle
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
#' @export
#'
gumbelEstThreshold <- function(x, lt, n_min, n_max) {
  UseMethod('gumbelEstThreshold')
}

#' @export
gumbelEstThreshold.declustered_series <- function(x, lt,
                                                  n_min,
                                                  n_max) {

  gumbelEstThreshold.default(x = x$declustered_series,
                             lt = lt,
                             n_min = n_min,
                             n_max = n_max)
}

#' @export
gumbelEstThreshold.default <- function(x, lt, n_min, n_max) {

  thresholds <- genThresholds(x, n_min, n_max)
  w_stats <- rep(NA, length(thresholds))

  pb <- progress::progress_bar$new(total = length(thresholds), clear = FALSE,
                                   complete = '*')
  pb$tick(0)

  for (i in 1:length(thresholds)) {

    y <- x[x > thresholds[i]]

    pot_fit <- gumbelMLE.default(x = y,
                                 lt = lt,
                                 thresh = thresholds[i],
                                 hessian_tf = FALSE)
    w_stats[i] <- gumbelWPlot.default(x = pot_fit$par,
                                      y = y,
                                      thresh = thresholds[i],
                                      tf_plot = FALSE,
                                      BW = FALSE,
                                      details = FALSE)

    pb$tick()
  }

  thresh <- thresholds[w_stats == min(w_stats)]

  value <- list(selected_threshold = thresh,
                lt = lt,
                y = x[x > thresh],
                checked_thresholds = thresholds,
                w_stats = w_stats)
  class(value) <- 'thresholded_series'
  value
}

#'
#' @title fullEstThreshold
#'
#' @description fullEstThreshold
#'
#' @details Select a threshold when using the full POT model
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
#' @export
#'
fullEstThreshold <- function (x, lt, n_min, n_max, n_starts) {
  UseMethod('fullEstThreshold')
}

#' @export
fullEstThreshold.declustered_series <- function(x, lt,
                                                n_min,
                                                n_max,
                                                n_starts) {

  fullEstThreshold.default(x = x$declustered_series,
                           lt = lt,
                           n_min = n_min,
                           n_max = n_max,
                           n_starts = n_starts)
}

#' @export
fullEstThreshold.default <- function(x, lt, n_min, n_max,
                                     n_starts) {

  thresholds <- genThresholds(x, n_min, n_max)
  w_stats <- rep(NA, length(thresholds))

  pb <- progress::progress_bar$new(total = length(thresholds), clear = FALSE,
                                   complete = '*')
  pb$tick(0)

  for (i in 1:length(thresholds)) {

    y <- x[x > thresholds[i]]

    pot_fit <- fullMLE.default(x = y,
                               lt = lt,
                               thresh = thresholds[i],
                               n_starts = n_starts,
                               hessian_tf = FALSE)
    w_stats[i] <- fullWPlot.default(x = pot_fit$par,
                                    y = y,
                                    thresh = thresholds[i],
                                    tf_plot = FALSE,
                                    BW = FALSE,
                                    details = FALSE)

    pb$tick()
  }

  thresh <- thresholds[w_stats == min(w_stats)]

  value <- list(selected_threshold = thresh,
                lt = lt,
                y = x[x > thresh],
                checked_thresholds = thresholds,
                w_stats = w_stats)
  class(value) <- 'thresholded_series'
  value
}
