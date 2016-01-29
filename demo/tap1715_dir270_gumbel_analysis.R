complete_series <- -jp1tap1715wind270$value

declustered_obs <- decluster(complete_series)

total_time <- 100 # seconds

thresholded_obs <- gumbelEstThreshold(x = declustered_obs,
                                      lt = total_time,
                                      n_min = 10,
                                      n_max = 100)

summary(thresholded_obs)

gumbel_pot_fit <- gumbelMLE(x = thresholded_obs,
                            hessian_tf = TRUE)

gumbelWPlot(x = gumbel_pot_fit,
            tf_plot = TRUE, BW = FALSE, details = FALSE)

gumbel_max_dist <- gumbelMaxDist(x = gumbel_pot_fit,
                                 lt_gen = 100,
                                 n_mc = 1000)
plot(gumbel_max_dist)

mean(gumbel_max_dist)

gumbel_max_dist_uncert <- gumbelMaxDistUncert(x = gumbel_pot_fit,
                                              lt_gen = 100, n_mc = 1000, n_boot = 1000)

summary(gumbel_max_dist_uncert)

plot(gumbel_max_dist_uncert)
plot(gumbel_max_dist, add = TRUE, border = 'red')

gumbelAnalysis(complete_series = complete_series,
               length_series = total_time,
               n_min = 10,
               n_max = 100,
               length_target_series = c(100, 167),
               wplot_filename = './demo/tap1715_dir270_gumbel_wplot',
               max_dist_hist_filename = './demo/tap1716_dir270_gumbel_max_dist_hist')
