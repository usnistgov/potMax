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
# plot.full_max_dist <- function(x,
#                                add_mean = TRUE,
#                                mean_col = 'red',
#                                binwidth = NULL,
#                                mean_size = 5,
#                                ...) {
#   if (is.null(binwidth)) {
#     binwidth <- (max(x$max_dist) - min(x$max_dist))/30
#   }
#   p <- ggplot2::ggplot(data = data.frame(x = x$max_dist),
#                        mapping = ggplot2::aes(x = x, y = ..density..)) +
#     ggplot2::geom_histogram(binwidth = binwidth) +
#     ggplot2::xlab('Peak Value') +
#     ggplot2::ylab('Density') +
#     ggplot2::ggtitle('Distribution of the Peak Value') +
#     ggplot2::geom_point(data = data.frame(x = mean(x), y = 0),
#                         mapping = ggplot2::aes(x = x, y = y),
#                         col = mean_col, size = mean_size) +
#     ggplot2::theme_classic()
#   p
# }

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
  if (add) {
    p <- ggplot2::last_plot() +
      ggplot2::geom_density(data = ggdata,
                            mapping = ggplot2::aes(x = x, group = dist_indices),
                            color = ggplot2::alpha('black', dens_line_alpha))
  } else {
    p <- ggplot2::ggplot(data = ggdata,
                         mapping = ggplot2::aes(x = x, group = dist_indices)) +
      ggplot2::geom_density(color = ggplot2::alpha('black', dens_line_alpha))
  }

  p <- p + ggplot2::xlab('Peak Value') +
      ggplot2::ylab('Density') +
      ggplot2::ggtitle('Distribution of the Peak Value') +
      ggplot2::theme_classic()

  if (add_int) {

    int_df <- data.frame(x = unlist(summary(x, probs = c(0.5 - int_prob/2, 0.5 + int_prob/2),
                                            suppress = TRUE)[2:3]),
                         y = rep(0, times = 2))
    p <- p + ggplot2::geom_line(data = int_df,
                                mapping = ggplot2::aes(x = x, y = y, group = NULL),
                                col = int_col, size = int_size)
  }

  p
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
# plot.full_max_dist_uncert <- function(x,
#                                       add_int = TRUE,
#                                       int_col = 'red',
#                                       int_size = 3,
#                                       int_prob = 0.8,
#                                       n_plot = 50,
#                                       add = FALSE,
#                                       dens_line_alpha = 0.25,
#                                         ...) {
#
#   dist_indices <- sample(x = 1:nrow(x$boot_samps),
#                          size = n_plot, replace = FALSE)
#   ggdata <- data.frame(x = as.numeric(t(x$boot_samps[dist_indices, ])),
#                        dist_indices = rep(dist_indices, each = ncol(x$boot_samps)))
#   if (add) {
#     p <- ggplot2::last_plot() +
#       ggplot2::geom_density(data = ggdata,
#                             mapping = ggplot2::aes(x = x, group = dist_indices),
#                             color = ggplot2::alpha('black', dens_line_alpha))
#   } else {
#     p <- ggplot2::ggplot(data = ggdata,
#                          mapping = ggplot2::aes(x = x, group = dist_indices)) +
#       ggplot2::geom_density(color = ggplot2::alpha('black', dens_line_alpha))
#   }
#
#   p <- p + ggplot2::xlab('Peak Value') +
#       ggplot2::ylab('Density') +
#       ggplot2::ggtitle('Distribution of the Peak Value') +
#       ggplot2::theme_classic()
#
#   if (add_int) {
#
#     int_df <- data.frame(x = unlist(summary(x, probs = c(0.5 - int_prob/2, 0.5 + int_prob/2),
#                                             suppress = TRUE)[2:3]),
#                          y = rep(0, times = 2))
#     p <- p + ggplot2::geom_line(data = int_df,
#                                 mapping = ggplot2::aes(x = x, y = y, group = NULL),
#                                 col = int_col, size = int_size)
#   }
#
#   p
# }