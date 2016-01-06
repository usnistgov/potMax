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
fullWPlot <- function (mu,
                       sigma,
                       k,
                       y,
                       thresh,
                       tf_plot,
                       BW,
                       details,
                       ...) {

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
  W <- sort(x=W)

  ## calculate the appropriate exp(1)
  ## quantiles
  quantiles <- ((1:n) - 0.375)/(n + 0.25)
  exp1_quantiles <- qexp(p=quantiles)

  if (tf_plot) {

    plot(x=exp1_quantiles, y=W,
         xlab="exp(1) quantiles",
         ylab="Ordered W-Statistics",
         main="W-Statistic Plot",
         bty = "l",
         ...)

    if (BW) {

      abline(a=0, b=1, col="black")
    } else {

      abline(a=0, b=1, col="red")
    }
  }

  if (details) {

    value <- abs((W - exp1_quantiles))
    data.frame(y, excesses, W1, W2, W3, W4, W, value)
  } else {

    max(abs((W - exp1_quantiles)))
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
gumbelWPlot <- function (mu,
                         sigma,
                         y,
                         thresh,
                         tf_plot,
                         BW,
                         details,
                         ...) {

  fullWPlot(mu = mu,
            sigma = sigma,
            k = 0,
            y = y,
            thresh = thresh,
            tf_plot = tf_plot,
            BW = BW,
            details = details,
            ...)
}

#' @title genThresholds
#'
#' @description genThresholds
#'
#' @details Generate the thresholds from which the estimated threshold will be
#'   chosen
#'
#' @param y_all
#'
#' @param n_min
#'
#' @param n_max
#'
genThresholds <- function (y_all, n_min, n_max) {

  y <- sort(y_all, decreasing = TRUE)
  thresholds <- NULL

  for (i in 1:(length(y) - 1)) {

    tmp_thresh <- mean(c(y[i], y[i + 1]))
    tmp_y <- y[y > tmp_thresh]
    if (length(tmp_y) >= n_min && length(tmp_y) <= n_max) {

      thresholds <- c(thresholds, tmp_thresh)
    }
  }

  thresholds
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
gumbelEstThreshold <- function (y_all, lt, n_min, n_max) {

  thresholds <- genThresholds(y_all, n_min, n_max)
  w_stats <- rep(NA, length(thresholds))

  for (i in 1:length(thresholds)) {

    y <- y_all[y_all > thresholds[i]]

    pot_fit <- gumbelMLE(y, lt, thresholds[i],
                         hessian = FALSE)
    w_stats[i] <- gumbelWPlot(pot_fit$par[1],
                              pot_fit$par[2],
                              y,
                              thresholds[i],
                              tf_plot = FALSE,
                              BW = FALSE,
                              details = FALSE)
  }

  thresh = thresholds[w_stats == min(w_stats)]

  list(selected_threshold = thresh,
       y = y_all[y_all > thresh],
       checked_thresholds = thresholds,
       w_stats = w_stats)
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
fullEstThreshold <- function (y_all, lt, n_min, n_max) {

  thresholds <- genThresholds(y_all, n_min, n_max)
  w_stats <- rep(NA, length(thresholds))

  for (i in 1:length(thresholds)) {

    y <- y_all[y_all > thresholds[i]]

    pot_fit <- fullMLE(y, lt, thresholds[i],
                       hessian = FALSE)
    w_stats[i] <- fullWPlot(pot_fit$par[1],
                            pot_fit$par[2],
                            pot_fit$par[3],
                            y,
                            thresholds[i],
                            tf_plot = FALSE,
                            BW = FALSE,
                            details = FALSE)
  }

  thresh = thresholds[w_stats == min(w_stats)]

  list(selected_threshold = thresh,
       y = y_all[y_all > thresh],
       checked_thresholds = thresholds,
       w_stats = w_stats)
}
