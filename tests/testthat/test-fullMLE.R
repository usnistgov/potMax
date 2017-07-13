library(potMax)
context('Full Maximum Likelihood Estimate')

test_that('fullMLE produces expected result', {
  set.seed(123456)
  a <- decluster(-scan('../../data/jp1tap1715wind270.csv',
                       skip = 1, quiet = TRUE))
  a1 <- a$declustered_series[a$declustered_series > 1]
  mle <- fullMLE(x = a1, hessian_tf = TRUE,
                 lt = 100, thresh = 1, n_starts = 20)
  expect_equal(round(mle$par, 2),
               c(1.77, 0.37, 0.19))
  expect_equal(round(mle$hessian/100, 0)*100,
               matrix(data = c(-30200,
                               65300,
                               -29900,
                               65300,
                               -148500,
                               69500,
                               -29900,
                               69500,
                               -33600),
                      nrow = 3))
  set.seed(as.integer(Sys.time()))
})