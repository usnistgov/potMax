#' @export
summary.thresholded_series <- function(object) {

  cat('Selected threshold: ',
      object$selected_threshold,
      '\n\n')
  cat('Number of thresholded observations: ',
      length(object$y),
      '\n\n')
  cat('Summary statistics:\n')
  print(summary(object$y))
}

tabSetUp <- function(se, quants, preds, probs) {

  value <- data.frame(se, t(quants), t(preds))
  colnames(value) <- c('SE Mean',
                       paste0((probs[2] - probs[1])*100, '% LB Mean'),
                       paste0((probs[2] - probs[1])*100, '% UB Mean'),
                       paste0((probs[2] - probs[1])*100, '% LB PI'),
                       paste0((probs[2] - probs[1])*100, '% UB PI'))
  value
}

#' @export
summary.gumbel_max_dist_uncert <- function(object, probs = c(0.1, 0.9),
                                           suppress = FALSE) {

  quants <- quantile(apply(object$boot_samps, 1, mean), probs = probs)
  preds <- quantile(as.vector(object$boot_samps), probs = probs)
  se <- sd(apply(object$boot_samps, 1, mean))
  value <- tabSetUp(se = se, quants = quants,
                    preds = preds, probs = probs)
  if (!suppress) {
    print(value)
  }
  invisible(value)
}

#' @export
summary.full_max_dist_uncert <- function(object, probs = c(0.1, 0.9),
                                         suppress = FALSE) {
  summary.gumbel_max_dist_uncert(object = object,
                                 probs = probs,
                                 suppress = suppress)
}
# summary.full_max_dist_uncert <- function(object, probs = c(0.1, 0.9),
#                                          suppress = FALSE) {
#
#   quants <- quantile(apply(object$boot_samps, 1, mean), probs = probs)
#   preds <- quantile(as.vector(object$boot_samps), probs = probs)
#   se <- sd(apply(object$boot_samps, 1, mean))
#   value <- data.frame(se, t(quants), t(preds))
#   colnames(value) <- c('SE Mean',
#                        paste0((probs[2] - probs[1])*100, '% LB Mean'),
#                        paste0((probs[2] - probs[1])*100, '% UB Mean'),
#                        paste0((probs[2] - probs[1])*100, '% LB PI'),
#                        paste0((probs[2] - probs[1])*100, '% UB PI'))
#   if (!suppress) {
#     print(value)
#   }
#   invisible(value)
# }

#' @export
summary.gumbel_max_dist_uncert_multi_thresh <- function(object, probs = c(0.1, 0.9),
                                                        suppress = FALSE) {

  quants <- quantile(unlist(lapply(object, function(x)mean(x$max_dist))),
                     probs = probs)
  preds <- quantile(unlist(lapply(object, function(x)x$max_dist)),
                    probs = probs)
  se <- sd(unlist(lapply(object, function(x)mean(x$max_dist))))
  value <- tabSetUp(se = se, quants = quants,
                    preds = preds, probs = probs)
  if (!suppress) {
    print(value)
  }
  invisible(value)
}

#' @export
summary.full_max_dist_uncert_multi_thresh <- function(object, probs = c(0.1, 0.9),
                                                      suppress = FALSE) {
  summary.gumbel_max_dist_uncert_multi_thresh(object = object,
                                              probs = probs,
                                              suppress = suppress)
}

#' @export
mean.gumbel_max_dist <- function(x) {
  mean(x$max_dist)
}

#' @export
mean.full_max_dist <- function(x) {
  mean.gumbel_max_dist(x)
}

#' @export
mean.gumbel_max_dist_multi_thresh <- function(x) {
  mean.gumbel_max_dist(x)
}

#' @export
mean.full_max_dist_multi_thresh <- function(x) {
  mean.gumbel_max_dist(x)
}

# #' @export
meanAndSd <- function(x){
  UseMethod('meanAndSd')
}

# #' @export
meanAndSd.gumbel_max_dist <- function(x) {
  c(mean(x$max_dist), sd(x$max_dist))
}

# #' @export
meanAndSd.full_max_dist <- function(x) {
  meanAndSd.gumbel_max_dist(x)
}

# #' @export
meanAndSd.gumbel_max_dist_multi_thresh <- function(x) {
  meanAndSd.gumbel_max_dist(x)
}

# #' @export
meanAndSd.full_max_dist_multi_thresh <- function(x) {
  meanAndSd.gumbel_max_dist(x)
}

#' @export
summary.full_pot_fit <- function(object, suppress = FALSE) {

  value <- data.frame(t(object$par))
  colnames(value) <- c('mu', 'sigma', 'k')
  if (!suppress) {
   print(value)
  }
  invisible(value)
}

#' @export
summary.gumbel_pot_fit <- function(object, suppress = FALSE) {

  value <- data.frame(t(object$par))
  colnames(value) <- c('mu', 'sigma')
  if (!suppress) {
    print(value)
  }
  invisible(value)
}

#' @export
summary.gumbel_multi_fit <- function(object, suppress = FALSE) {

  mu <- sapply(object$all_fits, function(x)x$par[1])
  sigma <- sapply(object$all_fits, function(x)x$par[2])
  thresh <- sapply(object$all_fits, function(x)x$thresh)
  weights <- object$weights
  value <- data.frame(mu = mu,
                      sigma = sigma,
                      thresh = thresh,
                      weight = weights)
  if (!suppress) {
    print(value)
  }
  invisible(value)
}

#' @export
summary.full_multi_fit <- function(object, suppress = FALSE) {

  mu <- sapply(object$all_fits, function(x)x$par[1])
  sigma <- sapply(object$all_fits, function(x)x$par[2])
  k <- sapply(object$all_fits, function(x)x$par[3])
  thresh <- sapply(object$all_fits, function(x)x$thresh)
  weights <- object$weights
  value <- data.frame(mu = mu,
                      sigma = sigma,
                      k = k,
                      thresh = thresh,
                      weight = weights)
  if (!suppress) {
    print(value)
  }
  invisible(value)
}

#' @export
summary.gumbel_N_year_uncert_multi_thresh <- function(object, conf_int_level = 0.8) {

  N_year_vals <- sapply(object, function(x){x$N_year_val})
  list(standard_error = sd(N_year_vals),
       conf_int = quantile(N_year_vals, probs = c((1 - conf_int_level)/2,
                                                  (0.5 + conf_int_level/2))))
}