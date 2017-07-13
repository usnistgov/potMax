library(potMax)
context('Gumbel Maximum Likelihood Estimate')

test_that('gumbelMLE produces expected result', {
  a <- decluster(-scan('../../data/jp1tap1715wind270.csv',
                       skip = 1, quiet = TRUE))
  a1 <- a$declustered_series[a$declustered_series > 1]
  mle <- gumbelMLE(x = a1, hessian_tf = TRUE,
                   lt = 100, thresh = 1)
  expect_equal(round(mle$par, 3),
               c(1.723, 0.272))
  expect_equal(round(mle$hessian, 0),
               matrix(data = c(-19457, 51829, 51829, -157518),
                      nrow = 2))
})