
all_dat <- jp1tap813wind315$value
tdat <- -all_dat[1:floor(0.33*length(all_dat))]
vdat <- -all_dat[ceiling(0.33*length(all_dat)):length(all_dat)]

dtdat <- decluster(tdat)

best_thresh <- gumbelEstThreshold(x = dtdat, 33, n_min = 10, n_max = 50)

threshPlot(best_thresh)

best_est <- gumbelMLE(x = best_thresh, hessian_tf = TRUE)

multi_est <- gumbelMultiFit(x = dtdat, lt = 33, n_min = 10, n_max = 50, weight_scale = 5)

threshPlot(multi_est)

best_dist <- gumbelMaxDist(x = best_est, lt_gen = 67, n_mc = 5000)

multi_dist <- gumbelMaxDist(x = multi_est, lt_gen = 67, n_mc = 5000)

multi_dist_uncert <- gumbelMaxDistUncert(x = multi_est,
                                         declust_obs = dtdat$declustered_series,
                                         lt_gen = 67, n_mc = 5000, n_boot = 100, progress_tf = T)

best_dist_uncert <- gumbelMaxDistUncert(x = best_est, lt_gen = 67, n_mc = 5000, n_boot = 100)


plot(best_dist)
plot(best_dist_uncert, add = TRUE)
plot(multi_dist)
plot(multi_dist_uncert, add = TRUE)

fbest_thresh <- fullEstThreshold(x = dtdat, lt = 33, n_min = 20, n_max = 200, n_starts = 2)
fbest_est <- fullMLE(x = fbest_thresh, hessian_tf = TRUE, n_starts = 1)

threshPlot(fbest_thresh)

fmulti_est <- fullMultiFit(x = dtdat, lt = 33, n_min = 20, n_max = 200,
                           weight_scale = 5, n_starts = 1)

threshPlot(fmulti_est)

fbest_dist <- fullMaxDist(x = fbest_est, lt_gen = 67, n_mc = 5000)

fmulti_dist <- fullMaxDist(x = fmulti_est, lt_gen = 67, n_mc = 5000)

fmulti_dist_uncert <- fullMaxDistUncert(x = fmulti_est,
                                        declust_obs = dtdat$declustered_series,
                                        lt_gen = 67, n_mc = 5000, n_boot = 100,
                                        n_starts = 1, progress_tf = T)

fbest_dist_uncert <- fullMaxDistUncert(x = fbest_est, lt_gen = 67, n_mc = 5000, n_boot = 100)

plot(fbest_dist)
plot(fbest_dist_uncert, add = TRUE)
plot(fmulti_dist)
plot(fmulti_dist_uncert, add = TRUE)
