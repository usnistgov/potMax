#' @export
plot.gumbel_max_dist <- function(x,
                                 add_mean = TRUE,
                                 mean_col = 'red',
                                 binwidth = NULL,
                                 mean_size = 5,
                                 ...) {

  if (is.null(binwidth)) {
    binwidth <- (max(x$max_dist) - min(x$max_dist))/30
  }
  p <- ggplot2::ggplot(data = data.frame(x = x$max_dist),
              mapping = ggplot2::aes(x = x, y = ..density..)) +
    ggplot2::geom_histogram(binwidth = binwidth) +
    ggplot2::xlab('Peak Value') +
    ggplot2::ylab('Density') +
    ggplot2::ggtitle('Distribution of the Peak Value') +
    ggplot2::geom_point(data = data.frame(x = mean(x), y = 0),
               mapping = ggplot2::aes(x = x, y = y),
               col = mean_col, size = mean_size) +
    ggplot2::theme_classic()
  p
}

#' @export
plot.gumbel_max_dist_multi_thresh <- function(x,
                                              add_mean = TRUE,
                                              mean_col = 'red',
                                              binwidth = NULL,
                                              mean_size = 5,
                                              ...) {

  plot.gumbel_max_dist(x = x,
                       add_mean = add_mean,
                       mean_col = mean_col,
                       binwidth = binwidth,
                       mean_size = mean_size,
                       ...)
}

#' @export
plot.full_max_dist <- function(x,
                               add_mean = TRUE,
                               mean_col = 'red',
                               binwidth = NULL,
                               mean_size = 5,
                               ...) {
  plot.gumbel_max_dist(x = x,
                       add_mean = add_mean,
                       mean_col = mean_col,
                       binwidth = binwidth,
                       mean_size = mean_size,
                       ...)
}

#' @export
plot.full_max_dist_multi_thresh <- function(x,
                                            add_mean = TRUE,
                                            mean_col = 'red',
                                            binwidth = NULL,
                                            mean_size = 5,
                                            ...) {
  plot.gumbel_max_dist(x = x,
                       add_mean = add_mean,
                       mean_col = mean_col,
                       binwidth = binwidth,
                       mean_size = mean_size,
                       ...)
}

createUncertPlot <- function(boot_df, int_df,
                             add_int,
                             int_col,
                             int_size,
                             int_prob,
                             n_plot,
                             add,
                             dens_line_alpha) {

  if (add) {
    p <- ggplot2::last_plot() +
      ggplot2::geom_density(data = boot_df,
                            mapping = ggplot2::aes(x = x, group = dist_indices),
                            color = ggplot2::alpha('black', dens_line_alpha))
  } else {
    p <- ggplot2::ggplot(data = boot_df,
                         mapping = ggplot2::aes(x = x, group = dist_indices)) +
      ggplot2::geom_density(color = ggplot2::alpha('black', dens_line_alpha))
  }
  p <- p + ggplot2::xlab('Peak Value') +
      ggplot2::ylab('Density') +
      ggplot2::ggtitle('Distribution of the Peak Value') +
      ggplot2::theme_classic()

  if (add_int) {

    p <- p + ggplot2::geom_line(data = int_df,
                                mapping = ggplot2::aes(x = x, y = y, group = NULL),
                                col = int_col, size = int_size)
  }

  p
}

#' @export
plot.gumbel_max_dist_uncert <- function(x,
                                        add_int = TRUE,
                                        int_col = 'red',
                                        int_size = 3,
                                        int_prob = 0.8,
                                        n_plot = 50,
                                        add = FALSE,
                                        dens_line_alpha = 0.25,
                                        ...) {

  dist_indices <- sample(x = 1:nrow(x$boot_samps),
                         size = n_plot, replace = FALSE)
  ggdata <- data.frame(x = as.numeric(t(x$boot_samps[dist_indices, ])),
                       dist_indices = rep(dist_indices, each = ncol(x$boot_samps)))

  int_df <- NULL
  if (add_int) {

    int_df <- data.frame(x = unlist(summary(x, probs = c(0.5 - int_prob/2, 0.5 + int_prob/2),
                                            suppress = TRUE)[2:3]),
                         y = rep(0, times = 2))
  }

  createUncertPlot(boot_df = ggdata,
                   int_df = int_df,
                   add_int = add_int,
                   int_col = int_col,
                   int_size = int_size,
                   int_prob = int_prob,
                   n_plot = n_plot,
                   add = add,
                   dens_line_alpha = dens_line_alpha)
}

#' @export
#'
plot.gumbel_max_dist_uncert_multi_thresh <- function(x,
                                                     add_int = TRUE,
                                                     int_col = 'red',
                                                     int_size = 3,
                                                     int_prob = 0.8,
                                                     n_plot = 50,
                                                     add = FALSE,
                                                     dens_line_alpha = 0.25,
                                                     ...) {

  dist_indices <- sample(x = 1:length(x),
                         size = n_plot, replace = FALSE)
  dist_indices_n <- unlist(lapply(x, function(x)sum(x$n_each)))
  dist_indices_n <- dist_indices_n[dist_indices]
  boot_df <- data.frame(x = unlist(lapply(x[dist_indices], function(x)x$max_dist)),
                        dist_indices = rep(dist_indices, dist_indices_n))

  int_df <- NULL
  if (add_int) {

    int_df <- data.frame(x = unlist(summary(x, probs = c(0.5 - int_prob/2, 0.5 + int_prob/2),
                                            suppress = TRUE)[2:3]),
                         y = rep(0, times = 2))
  }
  createUncertPlot(boot_df = boot_df,
                   int_df = int_df,
                   add_int = add_int,
                   int_col = int_col,
                   int_size = int_size,
                   int_prob = int_prob,
                   n_plot = n_plot,
                   add = add,
                   dens_line_alpha = dens_line_alpha)
}

#' @export
#'
plot.full_max_dist_uncert <- function(x,
                                      add_int = TRUE,
                                      int_col = 'red',
                                      int_size = 3,
                                      int_prob = 0.8,
                                      n_plot = 50,
                                      add = FALSE,
                                      dens_line_alpha = 0.25,
                                        ...) {
  plot.gumbel_max_dist_uncert(x = x,
                              add_int = add_int,
                              int_col = int_col,
                              int_size = int_size,
                              int_prob = int_prob,
                              n_plot = n_plot,
                              add = add,
                              dens_line_alpha = dens_line_alpha,
                              ...)
}

#' @export
plot.full_max_dist_uncert_multi_thresh <- function(x,
                                                   add_int = TRUE,
                                                   int_col = 'red',
                                                   int_size = 3,
                                                   int_prob = 0.8,
                                                   n_plot = 50,
                                                   add = FALSE,
                                                   dens_line_alpha = 0.25,
                                                   ...) {
  plot.gumbel_max_dist_uncert_multi_thresh(x = x,
                                           add_int = add_int,
                                           int_col = int_col,
                                           int_size = int_size,
                                           int_prob = int_prob,
                                           n_plot = n_plot,
                                           add = add,
                                           dens_line_alpha = dens_line_alpha,
                                           ...)
}

#'
#' @title Plot Fit Statistics versus Thresholds
#'
#' @description This function creates a plot of of the fit statistics versus the
#'   associated thresholds
#'
#' @details The fit statistics are the maximum vertical distances of the points
#'   to the \eqn{45^\circ} line in the W-plots described in \code{gumbelWPlot}
#'   and \code{fullWPlot}
#'
#' @param x An object of class \code{thresholded_series},
#'   \code{gumbel_multi_fit}, or \code{full_multi_fit}
#'
#' @return A \code{ggplot} object
#'
#' @examples
#'
#' \dontrun{
#'   complete_series <- -jp1tap1715wind270$value
#'
#'   declustered_obs <- decluster(complete_series)
#'
#'   thresholded_series <- gumbelEstThreshold(x = declustered_obs,
#'                                            lt = 100,
#'                                            n_min = 10,
#'                                            n_max = 100)
#'
#'   threshPlot(thresholded_series)
#' }
#'
#' @export
#'
threshPlot <- function(x, ...) {
  UseMethod('threshPlot')
}

#' @describeIn threshPlot
#'
#' @export
#'
threshPlot.thresholded_series <- function(x, ...) {
  threshPlot.default(w_stats = x$w_stats,
                     thresholds = x$checked_thresholds)
}

#' @describeIn threshPlot
#'
#' @export
#'
threshPlot.gumbel_multi_fit <- function(x, ...) {
  threshPlot.default(w_stats = x$w_stats,
                     thresholds = x$thresholds)
}

#' @describeIn threshPlot
#'
#' @export
#'
threshPlot.full_multi_fit <- function(x, ...) {
  threshPlot.gumbel_multi_fit(x = x)
}

# #' @export
threshPlot.default <- function(w_stats, thresholds) {

  mindf <- data.frame(w_stats = min(w_stats),
                      thresholds = thresholds[w_stats == min(w_stats)])
  ggdf <- data.frame(w_stats = w_stats, thresholds = thresholds)
  p <- ggplot2::ggplot(data = ggdf,
                       mapping = ggplot2::aes(x = thresholds, y = w_stats)) +
    ggplot2::geom_point() +
    ggplot2::geom_point(data = mindf, color = 'red', size = 3) +
    ggplot2::xlab('Threshold') +
    ggplot2::ylab('W-statistic') +
    ggplot2::theme_classic()
  p
}