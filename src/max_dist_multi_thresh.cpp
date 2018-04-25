#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

// [[Rcpp::export]]
NumericVector gumbelMaxDistMultiCpp(NumericVector mu,
                                    NumericVector sigma,
                                    NumericVector Lambda,
                                    NumericVector integration_constant,
                                    int n_mc,
                                    bool progress_tf) {

  NumericVector max_dist(n_mc);
  int i, k;
  double n;
  double tmp;
  Progress p(n_mc, progress_tf);

  for (k = 0; k < n_mc; k++) {

    n = R::rpois(Lambda[k]);
    while (n < 1.0) {
      n = R::rpois(Lambda[k]);
    }

    max_dist[k] = mu[k] - sigma[k]*(log(integration_constant[k]) + log(1 - R::runif(0, 1)));

    for (i = 1; i < n; i++) {

      tmp = mu[k] - sigma[k]*(log(integration_constant[k]) + log(1 - R::runif(0, 1)));
      if (tmp > max_dist[k]) {

        max_dist[k] = tmp;
      }
    }

    p.increment();
  }

  return max_dist;
}
