#' @title W plot for the Full Poisson Process Intensity Function
#'
#' @description Plot of W-statistics versus standard exponential quantiles for
#'   the full Poisson process intensity function described in \code{fullMLE}
#'
#' @details This is the W plot from Chapter 1 of "Extreme Values in Finance,
#'   Telecomunications, and the Environment." The chapter was authored by
#'   Richard Smith. The formula for \eqn{W} is
#'
#'   \deqn{W = \frac{1}{k}\log\Big\{1 + \frac{ky}{\sigma + k(u - \mu)}\Big\}}
#'
#'   where \eqn{y} is the excess of (difference from) the threshold \eqn{u}
#'
#' @param x An S3 object of type \code{full_pot_fit} or numeric vector of length
#'   3 where the first, second, and thrid components are \eqn{\mu}, \eqn{\sigma},
#'   and \eqn{k}, respectively.
#'
#' @param tf_plot (logical scalar) Create the plot if TRUE, else not.
#'
#' @param BW (logical scalar) The plot is created in black and white if TRUE,
#'   else not
#'
#' @param details (logical scalar) Should details of the calculation be returned
#'   as a \code{data.frame}
#'
#' @return If details is TRUE, a \code{data.frame} with the details of the
#'   calculation, else invisibly return the maximum vertical distance from the
#'   points of the plot to the \eqn{45^\circ} line
#'
#' @examples
#'
#' \dontrun{
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
#'
#' thresholded_obs <- fullEstThreshold(x = declustered_obs,
#'                                     lt = 100,
#'                                     n_min = 10,
#'                                     n_max = 100)
#'
#' full_pot_fit <- fullMLE(x = thresholded_obs,
#'                         hessian_tf = TRUE)
#'
#' fullWPlot(x = full_pot_fit, tf_plot = TRUE, BW = FALSE, details = FALSE)
#' }
#'
#' @export
#'
fullWPlot <- function (x, tf_plot, BW, details, ...) {
  UseMethod('fullWPlot')
}

#' @describeIn fullWPlot
#'
#' @export
#'
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

#' @describeIn fullWPlot
#'
#' @param y (numeric vector) The observations that exceed the threshold, NOT the
#'   excesses or differences, but the actual observations
#'
#' @export
#'
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

#' @title W plot for the Gumbel Poisson Process Intensity Function
#'
#' @description Plot of W-statistics versus standard exponential quantiles for
#'   the Gumbel like intensity function described in \code{gumbelMLE}
#'
#' @details This is the W plot from Chapter 1 of "Extreme Values in Finance,
#'   Telecomunications, and the Environment," but modified by taking the limit
#'   of the formula for \eqn{W} as the tail length parameter goes to zero
#'
#'   \deqn{W = \frac{y}{\sigma}}
#'
#'   where \eqn{y} is the excess of (difference from) the threshold.
#'
#' @param x An S3 object of type \code{gumbel_pot_fit} or numeric vector of length
#'   2 where the first and second components are \eqn{\mu}, \eqn{\sigma},
#'   respectively.
#'
#' @param tf_plot (logical scalar) Create the plot if TRUE, else not.
#'
#' @param BW (logical scalar) The plot is created in black and white if TRUE,
#'   else not
#'
#' @param details (logical scalar) Should details of the calculation be returned
#'   as a \code{data.frame}
#'
#' @return If details is TRUE, a \code{data.frame} with the details of the
#'   calculation, else invisibly return the maximum vertical distance from the
#'   points of the plot to the \eqn{45^\circ} line
#'
#' @examples
#'
#' \dontrun{
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
#'
#' thresholded_obs <- gumbelEstThreshold(x = declustered_obs,
#'                                      lt = 100,
#'                                      n_min = 10,
#'                                      n_max = 100)
#'
#' gumbel_pot_fit <- gumbelMLE(x = thresholded_obs,
#'                             hessian_tf = TRUE)
#'
#' gumbelWPlot(x = gumbel_pot_fit, tf_plot = TRUE, BW = FALSE, details = FALSE)
#' }
#'
#' @export
#'
gumbelWPlot <- function (x, tf_plot, BW, details, ...) {
  UseMethod('gumbelWPlot')
}

#' @describeIn gumbelWPlot
#'
#' @export
#'
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

#' @describeIn gumbelWPlot
#'
#' @param y (numeric vector) The observations that exceed the threshold, NOT the
#'   excesses or differences, but the actual observations
#'
#' @export
#'
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

  for (i in seq_along(thresholds)) {

    tmp <- y[y > thresholds[i]]
    if (length(unique(tmp)) < 2) {
      thresholds[i] <- NA
    } else {
      break
    }
  }
  thresholds <- thresholds[!is.na(thresholds)]

  if (length(thresholds) < 1) {
    tmp <- unique(y)
    if (length(tmp) > 2) {
      thresholds <- mean(tmp[2:3])
    } else if (length(tmp) == 2) {
      thresholds <- 0.99*tmp[2]
    } else {
      stop(paste0('potMax:::genThresholds: The data are constantly ', tmp))
    }
  }

  thresholds
}

#' @title Find the Optimal Threshold for the Gumbel Model
#'
#' @description Estimate the threshold to use for the 2D extremal Poisson
#'   process when the tail length parameter for that model is exactly zero.
#'
#' @details A sequence of candidate thresholds is generated, and the threshold
#'   that minimizes the maximum vertical distance of the points of the W plot
#'   (described in \code{gumbelWplot}) to the \eqn{45^\circ} line is selected.
#'   See the vingette for more details.
#'
#' @param x An S3 object of class \code{declustered_series} or a numeric
#'   vector.
#'
#' @param lt (numeric scalar) The length of the time series in units
#'   of time (seconds, minutes, hours, etc.).
#'
#' @param n_min The minimum number of thresholded observations to include
#'
#' @param n_max The maximum number of thresholded observations to include
#'
#' @param progress_tf Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{thresholded_series} with elements
#'   \code{$seleted_threshold}, \code{$lt}, \code{$y},
#'   \code{$checked_thresholds}, and \code{$w_stats}.  The element \code{$y} is
#'   a numeric vector containing the actual observations (NOT differences from
#'   the threshold) that exceed the selected threshold.  The element
#'   \code{$w_stats} contains the maximum vertical distance from points to the
#'   \eqn{45^\circ} line of the W plot for the corresponding
#'   \code{$checked_threshold}.
#'
#' @examples
#'
#' \dontrun{
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
#'
#' thresholded_obs <- gumbelEstThreshold(x = declustered_obs,
#'                                      lt = 100,
#'                                      n_min = 10,
#'                                      n_max = 100)
#' }
#'
#' @export
#'
gumbelEstThreshold <- function(x, lt, n_min, n_max, progress_tf = TRUE) {
  UseMethod('gumbelEstThreshold')
}

#' @describeIn gumbelEstThreshold
#'
#' @export
#'
gumbelEstThreshold.declustered_series <- function(x, lt,
                                                  n_min,
                                                  n_max,
                                                  progress_tf = TRUE) {

  gumbelEstThreshold.default(x = x$declustered_series,
                             lt = lt,
                             n_min = n_min,
                             n_max = n_max,
                             progress_tf = progress_tf)
}

#' @describeIn gumbelEstThreshold
#'
#' @export
#'
gumbelEstThreshold.default <- function(x, lt, n_min, n_max,
                                       progress_tf = TRUE) {

  thresholds <- genThresholds(x, n_min, n_max)
  w_stats <- rep(NA, length(thresholds))

  if (progress_tf) {
    pb <- progress::progress_bar$new(total = length(thresholds), clear = FALSE,
                                     format = '|:bar| :percent ~ :eta',
                                     complete = '+', incomplete = ' ',
                                     current = ' ',
                                     width = floor(0.6*getOption('width')))

    pb$tick(0)
  }

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

    if (progress_tf) {
      pb$tick()
    }
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


#' @title Find the Optimal Threshold for the Full Model
#'
#' @description Estimate the threshold to use for the 2D extremal Poisson
#'   process.
#'
#' @details A sequence of candidate thresholds is generated, and the threshold
#'   that minimizes the maximum vertical distance of the points of the W plot
#'   (described in \code{fullWplot}) to the \eqn{45^\circ} line is selected. See
#'   the vingette for more details.
#'
#' @param x An S3 object of class \code{declustered_series} or a numeric
#'   vector.
#'
#' @param lt (numeric scalar) The length of the time series in units
#'   of time (seconds, minutes, hours, etc.).
#'
#' @param n_min (numeric scalar) The minimum number of thresholded observations
#'   to include
#'
#' @param n_max (numeric scalar) The maximum number of thresholded observations
#'   to include
#'
#' @param n_starts (numeric scalar) An iterative algorithm is used to calculate
#'   the MLE, and the optimization algorithm is run \code{n_starts} times,
#'   feeding in different random starts each time.  There is no default because
#'   setting it too high results in extended run times, but it should be at
#'   least one (obviously).
#'
#' @param progress_tf (logical scalar) Display a progress bar if TRUE, else not.
#'
#' @return An S3 object of class \code{thresholded_series} with elements
#'   \code{$seleted_threshold}, \code{$lt}, \code{$y},
#'   \code{$checked_thresholds}, and \code{$w_stats}.  The element \code{$y} is
#'   a numeric vector containing the actual observations (NOT differences from
#'   the threshold) that exceed the selected threshold.  The element
#'   \code{$w_stats} contains the maximum vertical distance from points to the
#'   \eqn{45^\circ} line of the W plot for the corresponding
#'   \code{$checked_threshold}.
#'
#' @examples
#'
#' \dontrun{
#' complete_series <- -jp1tap1715wind270$value
#'
#' declustered_obs <- decluster(complete_series)
#'
#' thresholded_obs <- fullEstThreshold(x = declustered_obs,
#'                                     lt = 100,
#'                                     n_min = 10,
#'                                     n_max = 100,
#'                                     n_starts = 10)
#' }
#'
#' @export
#'
fullEstThreshold <- function(x, lt, n_min, n_max, n_starts, progress_tf = TRUE) {
  UseMethod('fullEstThreshold')
}

#' @describeIn fullEstThreshold
#'
#' @export
#'
fullEstThreshold.declustered_series <- function(x, lt,
                                                n_min,
                                                n_max,
                                                n_starts,
                                                progress_tf = TRUE) {

  fullEstThreshold.default(x = x$declustered_series,
                           lt = lt,
                           n_min = n_min,
                           n_max = n_max,
                           n_starts = n_starts,
                           progress_tf = progress_tf)
}

#' @describeIn fullEstThreshold
#'
#' @export
#'
fullEstThreshold.default <- function(x, lt, n_min, n_max,
                                     n_starts, progress_tf = TRUE) {

  thresholds <- genThresholds(x, n_min, n_max)
  w_stats <- rep(NA, length(thresholds))

  if (progress_tf) {
    pb <- progress::progress_bar$new(total = length(thresholds), clear = FALSE,
                                     format = '|:bar| :percent ~ :eta',
                                     complete = '+', incomplete = ' ',
                                     current = ' ',
                                     width = floor(0.6*getOption('width')))
    pb$tick(0)
  }

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

    if (progress_tf) {
      pb$tick()
    }
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
