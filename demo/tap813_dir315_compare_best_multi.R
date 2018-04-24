library(ggplot2)

all_dat <- jp1tap813wind315$value
tdat <- -all_dat[1:floor(0.33*length(all_dat))]
vdat <- -all_dat[ceiling(0.33*length(all_dat)):length(all_dat)]

dtdat <- decluster(tdat)

best_thresh <- gumbelEstThreshold(x = dtdat, 33, n_min = 10, n_max = 50)
best_est <- gumbelMLE(x = best_thresh, hessian_tf = TRUE)

multi_est <- gumbelMultiFit(x = dtdat, lt = 33, n_min = 10, n_max = 50, weight_scale = 5)

best_dist <- gumbelMaxDist(x = best_est, lt_gen = 67, n_mc = 5000)

multi_dist <- gumbelMaxDist(x = multi_est, lt_gen = 67, n_mc = 5000)

system.time(multi_dist_uncert <- gumbelMaxDistUncert(x = multi_est,
                                         declust_obs = dtdat$declustered_series,
                                         lt_gen = 67, n_mc = 5000, n_boot = 100, progress_tf = T))

system.time(best_dist_uncert <- gumbelMaxDistUncert(x = best_est, lt_gen = 67, n_mc = 5000, n_boot = 100))

ggdat <- tibble(y = c(best_dist$max_dist, multi_dist$max_dist),
                type = rep(c('best', 'multi'),
                           c(length(best_dist$max_dist), length(multi_dist$max_dist))))

ggplot(ggdat, aes(x = y, fill = type)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = mean(best_dist$max_dist)) +
  geom_vline(xintercept = mean(multi_dist$max_dist), linetype = 'dashed') +
  geom_vline(xintercept = max(vdat), linetype = 'dotted') +
  theme_classic()

