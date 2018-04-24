
#' @export
gumbelMaxDist.gumbel_multi_fit <- function(x,
                                           lt_gen, n_mc,
                                           progress_tf = TRUE) {

  mu <- sapply(x$all_fits, function(x) x$par[1])
  sigma <- sapply(x$all_fits, function(x)x$par[2])
  thresh <- sapply(x$all_fits, function(x)x$thresh)
  lt_gen <- lt_gen[1] # only scalars are allowed

  if (sum(sigma <= 0) > 0) {
    stop('must have all sigma > 0')
  }
  if (lt_gen <= 0) {
    stop('must have lt_gen > 0')
  }

  Lambda <- lt_gen*exp((-(thresh - mu))/sigma)
  const <- Lambda/lt_gen

  prob_zero_obs <- dpois(0, Lambda)
  if (sum(prob_zero_obs > 0.9) > 0) {

    warning(paste0('The probability of zero observations is more than 0.9 for ', sum(prob_zero_obs > 0.9), ' thresholds\n'), immediate. = TRUE)
  }

  n_each <- round(n_mc*x$weights)
  long_mu <- rep(mu, n_each)
  long_sigma <- rep(sigma, n_each)
  long_thresh <- rep(thresh, n_each)
  long_Lambda <- rep(Lambda, n_each)
  long_const <- rep(const, n_each)

  value <- list(mu = mu,
                sigma = sigma,
                thresh = thresh,
                lt_gen = lt_gen,
                n_each = n_each,
                max_dist = gumbelMaxDistMultiCpp(long_mu, long_sigma,
                                                 long_Lambda,
                                                 long_const,
                                                 sum(n_each),
                                                 progress_tf))
  class(value) <- 'gumbel_max_dist_multi_thresh'
  value
}

#' @export
gumbelMaxDistUncert.gumbel_multi_fit <- function(x,
                                                 declust_obs,
                                                 lt_gen,
                                                 n_mc,
                                                 n_boot,
                                                 progress_tf = TRUE) {

  boot_samps <- lapply(1:n_boot, function(i, x)sample(x,length(x),TRUE),
                       x = declust_obs)

  boot_multi_fits <- lapply(boot_samps, gumbelMultiFit.default,
                            lt = x$lt,
                            n_min = x$n_min,
                            n_max = x$n_max,
                            weight_scale = x$weight_scale)

  if (progress_tf) {
    boot_max_dist <- pbapply::pblapply(boot_multi_fits, gumbelMaxDist,
                                       lt_gen = lt_gen, n_mc = n_mc,
                                       progress_tf = FALSE)
  } else {
    boot_max_dist <- lapply(boot_multi_fits, gumbelMaxDist,
                            lt_gen = lt_gen, n_mc = n_mc,
                            progress_tf = FALSE)
  }
  class(boot_max_dist) <- 'gumbel_max_dist_uncert_multi_thresh'
}