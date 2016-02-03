complete_series <- -jp1tap1715wind270$value

declustered_obs <- decluster(complete_series)

total_time <- 100 # seconds

thresholded_obs <- fullEstThreshold(x = declustered_obs,
                                    lt = total_time,
                                    n_starts = 20,
                                    n_min = 10,
                                    n_max = 50)

summary(thresholded_obs)

full_pot_fit <- fullMLE(x = thresholded_obs,
                        n_starts = 100,
                        hessian_tf = TRUE)

fullWPlot(x = full_pot_fit,
          tf_plot = TRUE, BW = FALSE, details = FALSE)

full_max_dist <- fullMaxDist(x = full_pot_fit,
                             lt_gen = 100,
                             n_mc = 1000)

mean(full_max_dist)

full_max_dist_uncert <- fullMaxDistUncert(x = full_pot_fit,
                                          lt_gen = 167, n_mc = 1000, n_boot = 1000)

summary(full_max_dist_uncert)

plot(full_max_dist)
plot(full_max_dist_uncert, add = TRUE)

fullAnalysis(complete_series = complete_series,
             length_series = total_time,
             n_min = 10,
             n_max = 100,
             length_target_series = c(100, 167),
             wplot_filename = './demo/full_wplot',
             max_dist_hist_filename = './demo/full_max_dist')
