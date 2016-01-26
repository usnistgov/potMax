#' @export
plot.gumbel_max_dist <- function (x, ...) {
  hist(x$max_dist,
       xlab = 'Maximum Pressure',
       main = 'Distribution of Maximum Pressure',
       ...)
}

#' @export
plot.full_max_dist <- function (x, ...) {
  hist(x$max_dist,
       xlab = 'Maximum Pressure',
       main = 'Distribution of Maximum Pressure',
       ...)
}