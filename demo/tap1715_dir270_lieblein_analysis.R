complete_series <- -jp1tap1715wind270$value

n_parts <- 4:50

results <- NULL

for (i in seq_along(n_parts)) {
  tmp <- liebBLUE(x = complete_series,
                  n_parts = n_parts[i],
                  target_n_parts = c(n_parts[i], ceiling(1.67*n_parts[i])),
                  probs = c(0.5704, 0.78))
  results <- rbind(results, c(n_parts[i], tmp[1, 1], tmp[2, 1], tmp[1, 2], tmp[2, 2]))
}

colnames(results) <- c('n', 'mean; relative length = 1',
                       'p = 0.78; relative length = 1',
                       'mean; relative length = 1.67',
                       'p = 0.78; relative length = 1.67')
write.csv(x = results, file = '../blue_results_to_Dat_2-10-2016.csv', row.names = FALSE)

tmp <- gumbelAnalysis(complete_series = complete_series,
                      length_series = 100, n_min = 10,
                      n_max = 100, length_target_series = c(100, 167))

plot(n_parts, results[, 2], type = 'l', xlab = 'n', ylab = 'Mean of max')
abline(h = tmp[1, 1], col = 'red', lty = 'dashed')

plot(n_parts, results[, 4], type = 'l', xlab = 'n', ylab = 'Mean of max')
abline(h = tmp[2, 1], col = 'red', lty = 'dashed')
