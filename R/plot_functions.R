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
    ggplot2::annotate('text', Inf, Inf,
                      label = paste0(mean_col, ' point is the mean'),
                      hjust = 1, vjust = 1) +
    ggplot2::theme_classic()
  p
}

#' @export
plot.full_max_dist <- function(x,
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
    ggplot2::annotate('text', Inf, Inf,
                      label = paste0(mean_col, ' point is the mean'),
                      hjust = 1, vjust = 1) +
    ggplot2::theme_classic()
  p
}

#' @export
plot.gumbel_max_dist_uncert <- function (x,
                                         add_int = TRUE,
                                         int_col = 'red',
                                         n_plot = 50,
                                         add = FALSE, ...) {

  dist_indices <- sample(x = 1:dim(x$boot_samps)[1],
                         size = n_plot, replace = FALSE)
  densx <- NULL
  densy <- NULL
  for (i in 1:n_plot) {

    tmp <- density(x$boot_samps[dist_indices[i], ])
    densx <- rbind(densx, tmp$x)
    densy <- rbind(densy, tmp$y)
  }

  if (add) {

    add_col <- rgb(128, 128, 128, 100,
                   maxColorValue = 255)
    for (i in 1:n_plot) {

      lines(densx[i, ], densy[i, ],
            col = add_col)
    }

    if (add_int) {

      lines(summary(x, suppress = TRUE)[2:3],
            rep(0, times = 2),
            col = int_col, lwd = 2)
      legend('right',
             legend = c('Bootstrap Replicates',
                        '80% CI for the Mean'),
             bty = 'n', col = c(add_col, int_col),
             lty = 'solid')
    } else {

      legend('right', legend = 'Bootstrap Replicates',
             bty = 'n', col = add_col, lty = 'solid')
    }
  } else {

    plot(c(min(densx), max(densx)),
         c(min(densy), max(densy)),
         ylab = 'Density',
         xlab = 'Peak Value',
         main = 'Bootstrap Replicates of the Distribution of the Peak',
         bty = 'l',
         type = 'n',
         ...)
    for (i in 1:n_plot) {

      lines(densx[i, ], densy[i, ], ...)
    }

    if (add_int) {

      lines(summary(x, suppress = TRUE)[2:3],
            rep(0, times = 2),
            col = int_col, lwd = 2)
      legend('topright', legend = '80% CI for Mean',
             col = int_col, lty = 'solid', bty = 'n')
    }
  }
}

#' @export
plot.full_max_dist_uncert <- function (x,
                                       add_int = TRUE,
                                       int_col = 'red',
                                       n_plot = 50,
                                       add = FALSE, ...) {

  dist_indices <- sample(x = 1:dim(x$boot_samps)[1],
                         size = n_plot, replace = FALSE)
  densx <- NULL
  densy <- NULL
  for (i in 1:n_plot) {

    tmp <- density(x$boot_samps[dist_indices[i], ])
    densx <- rbind(densx, tmp$x)
    densy <- rbind(densy, tmp$y)
  }

  if (add) {

    add_col <- rgb(128, 128, 128, 100,
                   maxColorValue = 255)
    for (i in 1:n_plot) {

      lines(densx[i, ], densy[i, ],
            col = add_col)
    }

    if (add_int) {

      lines(summary(x, suppress = TRUE)[2:3],
            rep(0, times = 2),
            col = int_col, lwd = 2)
      legend('right',
             legend = c('Bootstrap Replicates',
                        '80% CI for the Mean'),
             bty = 'n', col = c(add_col, int_col),
             lty = 'solid')
    } else {

      legend('right', legend = 'Bootstrap Replicates',
             bty = 'n', col = add_col, lty = 'solid')
    }
  } else {

    plot(c(min(densx), max(densx)),
         c(min(densy), max(densy)),
         ylab = 'Density',
         main = 'Bootstrap Replicates of the Distribution of the Peak',
         bty = 'l',
         type = 'n',
         ...)
    for (i in 1:n_plot) {

      lines(densx[i, ], densy[i, ], ...)
    }

    if (add_int) {

      lines(summary(x, suppress = TRUE)[2:3],
            rep(0, times = 2),
            col = int_col, lwd = 2)
      legend('topright', legend = '80% CI for Mean',
             col = int_col, lty = 'solid', bty = 'n')
    }
  }
}