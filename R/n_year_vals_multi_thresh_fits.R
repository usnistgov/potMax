#' @export
gumbelNYear.gumbel_multi_fit <- function(x, N) {

  mu <- sapply(x$all_fits, function(x) x$par[1])
  sigma <- sapply(x$all_fits, function(x)x$par[2])
  thresh <- sapply(x$all_fits, function(x)x$thresh)
  N <- N[1] # only scalars are allowed

  if (sum(sigma <= 0) > 0) {
    stop('must have all sigma > 0')
  }
  if (N < 0) {
    stop('must have N >= 1')
  }

  N_year_val <- sigma*log(N) + mu
  N_year_val <- sum(N_year_val*x$weights)

  value <- list(mu = mu,
                sigma = sigma,
                thresh = thresh,
                N = N,
                N_year_val = N_year_val)
  class(value) <- 'gumbel_N_year_val_multi_thresh'
  value
}

#' @export
gumbelNYearUncert.gumbel_multi_fit <- function(x,
                                               declust_obs,
                                               N,
                                               n_boot,
                                               progress_tf = TRUE) {

  if (progress_tf) {
    pb <- progress::progress_bar$new(total = n_boot, clear = FALSE,
                                     format = '|:bar| :percent ~ :eta',
                                     complete = '+', incomplete = ' ',
                                     current = ' ',
                                     width = floor(0.6*getOption('width')))
    pb$tick(0)
  }

  boot_samps <- lapply(1:n_boot, function(i, x)sample(x,length(x),TRUE),
                       x = declust_obs)

  boot_N_year <- list()

  for (i in 1:n_boot) {

    boot_multi_fit <- gumbelMultiFit.default(x = boot_samps[[i]],
                                             lt = x$lt,
                                             n_min = x$n_min,
                                             n_max = x$n_max,
                                             weight_scale = x$weight_scale,
                                             progress_tf = FALSE)
    boot_N_year[[i]] <- gumbelNYear(x = boot_multi_fit,
                                    N = N)
    if (progress_tf) {
      pb$tick()
    }
  }

  class(boot_N_year) <- 'gumbel_N_year_uncert_multi_thresh'
  boot_N_year
}
