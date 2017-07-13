library(potMax)
context('Generate Gumbel max dist, full max dist, Gumbel uncert and full uncert')

test_that('gumbelMaxDist produces expected result', {
  set.seed(123456)
  a <- decluster(-scan('../../data/jp1tap1715wind270.csv',
                       skip = 1, quiet = TRUE))
  a1 <- a$declustered_series[a$declustered_series > 1]
  mle <- gumbelMLE(x = a1, hessian_tf = TRUE,
                   lt = 100, thresh = 1)
  max_dist <- gumbelMaxDist(x = mle, lt_gen = 100, n_mc = 10,
                            progress_tf = FALSE)
  expect_equal(round(max_dist$max_dist, 5),
               c(2.82525,
                 3.37032,
                 3.52029,
                 3.24737,
                 3.17906,
                 3.00799,
                 3.08324,
                 2.95755,
                 3.34291,
                 3.02834))
  set.seed(as.integer(Sys.time()))
})

test_that('gumbelMaxDistUncert produces expected result', {
  set.seed(123456)
  a <- decluster(-scan('../../data/jp1tap1715wind270.csv',
                       skip = 1, quiet = TRUE))
  a1 <- a$declustered_series[a$declustered_series > 1]
  mle <- gumbelMLE(x = a1, hessian_tf = TRUE,
                   lt = 100, thresh = 1)
  max_dist_uncert <- gumbelMaxDistUncert(x = mle, lt_gen = 100,
                                         n_mc = 5, n_boot = 2,
                                         progress_tf = FALSE)
  expect_equal(round(max_dist_uncert$boot_samps, 5),
               matrix(data = c(2.79256,
                               3.32786,
                               3.47515,
                               3.20711,
                               3.14003,
                               3.02040,
                               3.09611,
                               2.96965,
                               3.35739,
                               3.04088),
                      nrow = 2, byrow = TRUE))
  set.seed(as.integer(Sys.time()))
})

test_that('fullMaxDist produces expected result', {
  set.seed(123456)
  a <- decluster(-scan('../../data/jp1tap1715wind270.csv',
                       skip = 1, quiet = TRUE))
  a1 <- a$declustered_series[a$declustered_series > 1]
  mle <- fullMLE(x = a1, hessian_tf = TRUE,
                 lt = 100, thresh = 1, n_starts = 20)
  max_dist <- fullMaxDist(x = mle, lt_gen = 100, n_mc = 10,
                          progress_tf = FALSE)
  expect_equal(round(max_dist$max_dist, 5),
               c(4.053040,
                 6.04723,
                 6.74633,
                 4.84280,
                 5.25766,
                 4.63719,
                 4.90081,
                 4.46827,
                 5.92735,
                 4.70710))
  set.seed(as.integer(Sys.time()))
})

test_that('fullMaxDistUncert produces expected result', {
  set.seed(123456)
  a <- decluster(-scan('../../data/jp1tap1715wind270.csv',
                       skip = 1, quiet = TRUE))
  a1 <- a$declustered_series[a$declustered_series > 1]
  mle <- fullMLE(x = a1, hessian_tf = TRUE,
                 lt = 100, thresh = 1, n_starts = 20)
  max_dist_uncert <- fullMaxDistUncert(x = mle, lt_gen = 100,
                                       n_mc = 5, n_boot = 2,
                                       progress_tf = FALSE)
  expect_equal(round(max_dist_uncert$boot_samps, 5),
               matrix(data = c(4.02788,
                               5.96266,
                               6.63607,
                               4.79687,
                               5.19930,
                               5.71303,
                               6.13523,
                               5.44714,
                               7.85907,
                               5.82413),
                      nrow = 2, byrow = TRUE))
  set.seed(as.integer(Sys.time()))
})