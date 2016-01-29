#' @export
plot.gumbel_max_dist <- function (x, ...) {
  hist(x$max_dist,
       xlab = 'Maximum Pressure',
       main = 'Distribution of Maximum Pressure',
       freq = FALSE,
       ...)
}

#' @export
plot.full_max_dist <- function (x, ...) {
  hist(x$max_dist,
       xlab = 'Maximum Pressure',
       main = 'Distribution of Maximum Pressure',
       freq = FALSE,
       ...)
}

#' @export
plot.gumbel_max_dist_uncert <- function (x, n_plot = 50,
                                         add = FALSE, ...) {

  dist_indices <- sample(x = 1:dim(x$boot_samps)[1],
                         size = n_plot, replace = FALSE)
  densx <- NULL
  densy <- NULL
  for (i in 1:n_plot) {

    tmp <- density(x$boot_samps[dist_indices[i], ], ...)
    densx <- rbind(densx, tmp$x)
    densy <- rbind(densy, tmp$y)
  }

  if (add) {

    for (i in 1:n_plot) {

      lines(densx[i, ], densy[i, ], ...)
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
  }
}

#' @export
plot.full_max_dist_uncert <- function (x, n_plot = 50,
                                       add = FALSE, ...) {

  dist_indices <- sample(x = 1:dim(x$boot_samps)[1],
                         size = n_plot, replace = FALSE)
  densx <- NULL
  densy <- NULL
  for (i in 1:n_plot) {

    tmp <- density(x$boot_samps[dist_indices[i], ], ...)
    densx <- rbind(densx, tmp$x)
    densy <- rbind(densy, tmp$y)
  }

  if (add) {

    lines(densx[i, ], densy[i, ], ...)
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
  }
}