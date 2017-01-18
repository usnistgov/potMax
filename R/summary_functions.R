#' @export
summary.thresholded_series <- function (object) {

  cat('Selected threshold: ',
      object$selected_threshold,
      '\n\n')
  cat('Number of thresholded observations: ',
      length(object$y),
      '\n\n')
  cat('Summary statistics:\n')
  print(summary(object$y))
}

#' @export
summary.gumbel_max_dist_uncert <- function(object, probs = c(0.1, 0.9),
                                           suppress = FALSE) {

  quants <- quantile(apply(object$boot_samps, 1, mean), probs = probs)
  preds <- quantile(as.vector(object$boot_samps), probs = probs)
  se <- sd(apply(object$boot_samps, 1, mean))
  value <- data.frame(se, t(quants), t(preds))
  colnames(value) <- c('SE Mean',
                       paste0((probs[2] - probs[1])*100, '% LB Mean'),
                       paste0((probs[2] - probs[1])*100, '% UB Mean'),
                       paste0((probs[2] - probs[1])*100, '% LB PI'),
                       paste0((probs[2] - probs[1])*100, '% UB PI'))
  if (!suppress) {
    print(value)
  }
  invisible(value)
}

#' @export
summary.full_max_dist_uncert <- function(object, probs = c(0.1, 0.9),
                                         suppress = FALSE) {

  quants <- quantile(apply(object$boot_samps, 1, mean), probs = probs)
  preds <- quantile(as.vector(object$boot_samps), probs = probs)
  se <- sd(apply(object$boot_samps, 1, mean))
  value <- data.frame(se, t(quants), t(preds))
  colnames(value) <- c('SE Mean',
                       paste0((probs[2] - probs[1])*100, '% LB Mean'),
                       paste0((probs[2] - probs[1])*100, '% UB Mean'),
                       paste0((probs[2] - probs[1])*100, '% LB PI'),
                       paste0((probs[2] - probs[1])*100, '% UB PI'))
  if (!suppress) {
    print(value)
  }
  invisible(value)
}

#' @export
mean.gumbel_max_dist <- function (x) {
  mean(x$max_dist)
}

#' @export
mean.full_max_dist <- function (x) {
  mean(x$max_dist)
}

#' @export
summary.full_pot_fit <- function (x) {

  value <- data.frame(t(x$par))
  colnames(value) <- c('mu', 'sigma', 'k')
  print(value)
  invisible(value)
}

#' @export
summary.gumbel_pot_fit <- function (x) {

  value <- data.frame(t(x$par))
  colnames(value) <- c('mu', 'sigma')
  print(value)
  invisible(value)
}