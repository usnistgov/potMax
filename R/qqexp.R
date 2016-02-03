# #'
# #' @title qqexp
# #'
# #' @description qqexp
# #'
# #' @details Generate a Q-Q plot for the exponential distribution
# #'
# #' @param y The data for which to generate the Q-Q plot
# #'
# #' @export
# #'
qqexp <- function (y, ...) {

  plot_y <- sort(y)/mean(y)
  plot_x_probs <- ((1:length(y)) - 0.375)/(length(y) + 0.25)
  plot_x <- qexp(p = plot_x_probs)
  plot(plot_x, plot_y,
       xlab = 'Theoretical Quantiles',
       ylab = 'Sample Quantiles',
       ...)
  abline(a = 0, b = 1)

  invisible(list(x = plot_x, y = plot_y))
}
