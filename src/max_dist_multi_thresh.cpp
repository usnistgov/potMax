#include <Rcpp.h>
using namespace Rcpp;

#include <RProgress.h>

// [[Rcpp::export]]
NumericVector gumbelMaxDistMultiCpp(NumericVector mu,
                                    NumericVector sigma,
                                    NumericVector Lambda,
                                    NumericVector integration_constant,
                                    int n_mc,
                                    bool progress_tf) {

  NumericVector max_dist(n_mc);
  int i, indx;
  double n;
  double tmp;
  RProgress::RProgress p("|:bar| :percent ~ :eta",
                         n_mc,
                         0.6*Rf_GetOptionWidth(),
                         '+',
                         ' ',
                         false,
                         0.2);

  for (indx = 0; indx < n_mc; indx++) {

    n = R::rpois(Lambda[indx]);
    while (n < 1.0) {
      n = R::rpois(Lambda[indx]);
    }

    max_dist[indx] = mu[indx] - sigma[indx]*(log(integration_constant[indx]) + log(1 - R::runif(0, 1)));

    for (i = 1; i < n; i++) {

      tmp = mu[indx] - sigma[indx]*(log(integration_constant[indx]) + log(1 - R::runif(0, 1)));
      if (tmp > max_dist[indx]) {

        max_dist[indx] = tmp;
      }
    }
    if (progress_tf) {
      p.tick();
    }
  }

  return max_dist;
}

// [[Rcpp::export]]
NumericVector fullMaxDistMultiCpp(NumericVector mu,
                                  NumericVector sigma,
                                  NumericVector k,
                                  NumericVector Lambda,
                                  NumericVector integration_constant,
                                  int n_mc,
                                  bool progress_tf) {

  NumericVector max_dist(n_mc);
  int i, indx;
  double n;
  double tmp;
  RProgress::RProgress p("|:bar| :percent ~ :eta",
                         n_mc,
                         0.6*Rf_GetOptionWidth(),
                         '+',
                         ' ',
                         false,
                         0.2);

  for (indx = 0; indx < n_mc; indx++) {

    n = R::rpois(Lambda[indx]);
    while (n < 1.0) {
      n = R::rpois(Lambda[indx]);
    }

    max_dist[indx] = (sigma[indx]/k[indx])*(pow(integration_constant[indx]*(1 - R::runif(0, 1)), -k[indx]) - 1) + mu[indx];

    for (i = 1; i < n; i++) {

      tmp = (sigma[indx]/k[indx])*(pow(integration_constant[indx]*(1 - R::runif(0, 1)), -k[indx]) - 1) + mu[indx];
      if (tmp > max_dist[indx]) {

        max_dist[indx] = tmp;
      }
    }
    if (progress_tf) {
      p.tick();
    }
  }

  return max_dist;
}
