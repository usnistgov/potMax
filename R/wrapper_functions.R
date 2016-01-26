#' @title gumbelAnalysis
#'
#' @description gumbelAnalysis
#'
#' @details gumbelAnalysis
#'
#' @param complete_series (numeric vector) The time series from the wind tunnel
#'   to be analyzed
#'
#' @param length_series (numeric scalar) The length of the time series in units
#'   of time (seconds, minutes, hours, etc.) from the wind tunnel
#'
#' @param n_min (numeric scalar) The minimum number of observations to include
#'   in the search for the optimal threshold
#'
#' @param m_max (numeric scalar) The maximum number of observations to include
#'   in the search for the optimal threshold
#'
#' @param length_target_series (numeric vector) The distribution of the peak
#'   depends on the length of observation, which could be different from the
#'   length of the original series.  Specify one or more lengths in units of
#'   time (seconds, minutes, hours, etc.) for the length of series for which the
#'   distribution of the peak is sought.  The units of time should be the same
#'   as the argument \code{length_series}.
#'
#' @param wplot_filename (NULL or string) If NULL the plot is not created; else,
#'   it is saved to "wplot_filename.pdf"
#'
#' @param BW (logical) Create the Wplot in black and white or color
#'
#' @param n_mc (numeric scalar) The number of Monte Carlo samples from which to
#'   construct the distribution of the peak
#'
#' @param max_dist_hist_filename (NULL or string) If NULL the plot is not
#'   created; else, it is saved to
#'   "max_dist_hist_filename_length_target_series.pdf"
#'
#' @param n_boot (numeric scalar) The number of bootstrap samples for the
#'   uncertainty evaluation
#'
#' @param ci_level (numeric scalar) The level of the confidence interval for the
#'   mean of the peak distribution.  Defaults to 0.8.
#'
#' @export
#'
gumbelAnalysis <- function (complete_series,
                            length_series,
                            n_min,
                            n_max,
                            length_target_series,
                            wplot_filename = NULL,
                            BW = NULL,
                            n_mc = 1000,
                            max_dist_hist_filename = NULL,
                            n_boot = 1000,
                            ci_level = 0.8) {

  cat('Start: declustering the series\n')
  declustered_obs <- decluster(complete_series)
  cat('Done: declustering te series \n\n')

  cat('Start: search for the optimal threshold\n')
  thresholded_obs <- gumbelEstThreshold(x = declustered_obs,
                                        lt = length_series,
                                        n_min = n_min,
                                        n_max = n_max)
  cat('Done: search for the optimal threshold\n\n')

  cat('Start: fitting Gumbel POT model to the optimal threshold\n')
  gumbel_pot_fit <- gumbelMLE(x = thresholded_obs,
                            hessian_tf = TRUE)
  cat('Done: fitting Gumbel POT model to the optimal threshold\n\n')

  if (!is.null(wplot_filename)) {

    cat(paste0('Start: saving wplot to file ', wplot_filename, '\n'))
    pdf(file = paste0(wplot_filename, '.pdf'))
    if (!is.null(BW)) {

      gumbelWPlot(x = gumbel_pot_fit,
                  tf_plot = TRUE, BW = BW, details = FALSE)
    } else {

      gumbelWPlot(x = gumbel_pot_fit,
                  tf_plot = TRUE, BW = TRUE, details = FALSE)
    }
    dev.off()
    cat(paste0('Done: saving wplot to file ', wplot_filename, '\n\n'))
  }

  cat('Start: estimating the distribution of the peak\n')
  gumbel_max_dist <- list(rep(NA, length(length_target_series)))
  for (i in seq_along(length_target_series)) {

    gumbel_max_dist[[i]] <- gumbelMaxDist(x = gumbel_pot_fit,
                                          lt_gen = length_target_series[i],
                                          n_mc = n_mc)
  }
  cat('Done: estimating the distribution of the peak\n\n')

  if (!is.null(max_dist_hist_filename)) {

    cat('Start: plotting the distribution of the peak\n')
    for (i in seq_along(length_target_series)) {

      pdf(paste0(max_dist_hist_filename, '_',
                 length_target_series[i], '.pdf'))
      plot(gumbel_max_dist[[i]])
      dev.off()
    }
    cat('Done: plotting the distribution of the peak\n\n')
  }

  cat('Start: estimating the uncertainty in the distribution of the peak\n')
  gumbel_max_dist_uncert <- list(rep(NA, length(length_target_series)))
  for (i in seq_along(length_target_series)) {

    gumbel_max_dist_uncert[[i]] <- gumbelMaxDistUncert(x = gumbel_pot_fit,
                                                       lt_gen = length_target_series[i],
                                                       n_mc = n_mc,
                                                       n_boot = n_boot)
  }
  cat('Done: estimating the uncertainty in the distribution of the peak\n\n')

  value <- matrix(nrow = length(length_target_series), ncol = 4)
  ci_probs <- c((1 - ci_level)/2, 1 - (1 - ci_level)/2)
  for (i in seq_along(length_target_series)) {

    value[i, 1] <- mean(gumbel_max_dist[[i]])
    value[i, 2:4] <- unlist(summary(gumbel_max_dist_uncert[[i]],
                                    probs = ci_probs,
                                    suppress = TRUE))
  }
  colnames(value) <- c('Mean', 'SE',
                       paste0(100*ci_level, '% LB'),
                       paste0(100*ci_level, '% UB'))
  print(value)
  invisible(value)
}

#' @title fullAnalysis
#'
#' @description fullAnalysis
#'
#' @details fullAnalysis
#'
#' @param complete_series (numeric vector) The time series from the wind tunnel
#'   to be analyzed
#'
#' @param length_series (numeric scalar) The length of the time series in units
#'   of time (seconds, minutes, hours, etc.) from the wind tunnel
#'
#' @param n_min (numeric scalar) The minimum number of observations to include
#'   in the search for the optimal threshold
#'
#' @param m_max (numeric scalar) The maximum number of observations to include
#'   in the search for the optimal threshold
#'
#' @param length_target_series (numeric vector) The distribution of the peak
#'   depends on the length of observation, which could be different from the
#'   length of the original series.  Specify one or more lengths in units of
#'   time (seconds, minutes, hours, etc.) for the length of series for which the
#'   distribution of the peak is sought.  The units of time should be the same
#'   as the argument \code{length_series}.
#'
#' @param n_starts (numeric scalar) The number of attempts to make at optimizing
#'   the likelihood when estimating the model parameters.  Setting this to a
#'   large value will make the algorithm, especially choosing the optimal
#'   threshold, take much longer.  Defaults to 20.
#'
#' @param wplot_filename (NULL or string) If NULL the plot is not created; else,
#'   it is saved to "wplot_filename.pdf"
#'
#' @param BW (logical) Create the Wplot in black and white or color
#'
#' @param n_mc (numeric scalar) The number of Monte Carlo samples from which to
#'   construct the distribution of the peak
#'
#' @param max_dist_hist_filename (NULL or string) If NULL the plot is not
#'   created; else, it is saved to
#'   "max_dist_hist_filename_length_target_series.pdf"
#'
#' @param n_boot (numeric scalar) The number of bootstrap samples for the
#'   uncertainty evaluation
#'
#' @param ci_level (numeric scalar) The level of the confidence interval for the
#'   mean of the peak distribution.  Defaults to 0.8.
#'
#' @export
#'
fullAnalysis <- function (complete_series,
                          length_series,
                          n_min,
                          n_max,
                          length_target_series,
                          n_starts = 20,
                          wplot_filename = NULL,
                          BW = NULL,
                          n_mc = 1000,
                          max_dist_hist_filename = NULL,
                          n_boot = 1000,
                          ci_level = 0.8) {

  cat('Start: declustering the series\n')
  declustered_obs <- decluster(complete_series)
  cat('Done: declustering te series \n\n')

  cat('Start: search for the optimal threshold\n')
  thresholded_obs <- fullEstThreshold(x = declustered_obs,
                                      lt = length_series,
                                      n_min = n_min,
                                      n_max = n_max,
                                      n_starts = n_starts)
  cat('Done: search for the optimal threshold\n\n')

  cat('Start: fitting full POT model to the optimal threshold\n')
  full_pot_fit <- fullMLE(x = thresholded_obs,
                          hessian_tf = TRUE,
                          n_starts = n_starts)
  cat('Done: fitting full POT model to the optimal threshold\n\n')

  if (!is.null(wplot_filename)) {

    cat(paste0('Start: saving wplot to file ', wplot_filename, '\n'))
    pdf(file = paste0(wplot_filename, '.pdf'))
    if (!is.null(BW)) {

      fullWPlot(x = full_pot_fit,
               tf_plot = TRUE, BW = BW, details = FALSE)
    } else {

      fullWPlot(x = full_pot_fit,
                tf_plot = TRUE, BW = TRUE, details = FALSE)
    }
    dev.off()
    cat(paste0('Done: saving wplot to file ', wplot_filename, '\n\n'))
  }

  cat('Start: estimating the distribution of the peak\n')
  full_max_dist <- list(rep(NA, length(length_target_series)))
  for (i in seq_along(length_target_series)) {

    full_max_dist[[i]] <- fullMaxDist(x = full_pot_fit,
                                      lt_gen = length_target_series[i],
                                      n_mc = n_mc)
  }
  cat('Done: estimating the distribution of the peak\n\n')

  if (!is.null(max_dist_hist_filename)) {

    cat('Start: plotting the distribution of the peak\n')
    for (i in seq_along(length_target_series)) {

      pdf(paste0(max_dist_hist_filename, '_',
                 length_target_series[i], '.pdf'))
      plot(full_max_dist[[i]])
      dev.off()
    }
    cat('Done: plotting the distribution of the peak\n\n')
  }

  cat('Start: estimating the uncertainty in the distribution of the peak\n')
  full_max_dist_uncert <- list(rep(NA, length(length_target_series)))
  for (i in seq_along(length_target_series)) {

    full_max_dist_uncert[[i]] <- fullMaxDistUncert(x = full_pot_fit,
                                                   lt_gen = length_target_series[i],
                                                   n_mc = n_mc,
                                                   n_boot = n_boot)
  }
  cat('Done: estimating the uncertainty in the distribution of the peak\n\n')

  value <- matrix(nrow = length(length_target_series), ncol = 4)
  ci_probs <- c((1 - ci_level)/2, 1 - (1 - ci_level)/2)
  for (i in seq_along(length_target_series)) {

    value[i, 1] <- mean(full_max_dist[[i]])
    value[i, 2:4] <- unlist(summary(full_max_dist_uncert[[i]],
                                    probs = ci_probs,
                                    suppress = TRUE))
  }
  colnames(value) <- c('Mean', 'SE',
                       paste0(100*ci_level, '% LB'),
                       paste0(100*ci_level, '% UB'))
  print(value)
  invisible(value)
}