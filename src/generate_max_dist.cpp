#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

// [[Rcpp::export]]
NumericVector gumbelMaxDistCpp(double mu, double sigma,
                               double Lambda,
                               double integration_constant,
                               int n_mc, bool progress_tf) {

  NumericVector max_dist(n_mc);
  int indx, i;
  double n;
  double tmp;
  Progress p(n_mc, progress_tf);

  for (indx = 0; indx < n_mc; indx++) {

    n = R::rpois(Lambda);
    while (n < 1.0) {
      n = R::rpois(Lambda);
    }

    max_dist[indx] = mu - sigma*(log(integration_constant) + log(1 - R::runif(0, 1)));

    for (i = 1; i < n; i++) {

      tmp = mu - sigma*(log(integration_constant) + log(1 - R::runif(0, 1)));
      if (tmp > max_dist[indx]) {

        max_dist[indx] = tmp;
      }
    }

    p.increment();
  }

  return max_dist;
}

// [[Rcpp::export]]
NumericMatrix gumbelMaxDistUncertCpp(NumericVector mu,
                                     NumericVector sigma,
                                     NumericVector Lambda,
                                     NumericVector integration_constant,
                                     int n_mc,
                                     int n_boot,
                                     bool progress_tf) {

  int indx, i, j;
  NumericMatrix max_dist_uncert(n_boot, n_mc);
  double n;
  double tmp;
  Progress p(n_boot, progress_tf);

  for (j = 0; j < n_boot; j++) {

    for (indx = 0; indx < n_mc; indx++) {

      n = R::rpois(Lambda[j]);
      while (n < 1.0) {
        n = R::rpois(Lambda[j]);
      }

      max_dist_uncert(j, indx) = mu[j] - sigma[j]*(log(integration_constant[j]) + log(1 - R::runif(0, 1)));

      for (i = 1; i < n; i++) {

        tmp = mu[j] - sigma[j]*(log(integration_constant[j]) + log(1 - R::runif(0, 1)));
        if (tmp > max_dist_uncert(j, indx)) {

          max_dist_uncert(j, indx) = tmp;
        }
      }
    }
    p.increment();
  }

  return max_dist_uncert;
}

// [[Rcpp::export]]
NumericVector fullMaxDistCpp(double mu, double sigma, double k,
                             double Lambda,
                             double integration_constant,
                             int n_mc,
                             bool progress_tf) {

  NumericVector max_dist(n_mc);
  int indx, i;
  double n;
  double tmp;
  Progress p(n_mc, progress_tf);

  for (indx = 0; indx < n_mc; indx++) {

    n = R::rpois(Lambda);
    while (n < 1.0){
      n = R::rpois(Lambda);
    }

    max_dist[indx] = (sigma/k)*(pow(integration_constant*(1 - R::runif(0, 1)), -k) - 1) + mu;

    for (i = 1; i < n; i++) {

      tmp = (sigma/k)*(pow(integration_constant*(1 - R::runif(0, 1)), -k) - 1) + mu;
      if (tmp > max_dist[indx]) {

        max_dist[indx] = tmp;
      }
    }
    p.increment();
  }

  return max_dist;
}

// [[Rcpp::export]]
NumericMatrix fullMaxDistUncertCpp(NumericVector mu,
                                   NumericVector sigma,
                                   NumericVector k,
                                   NumericVector Lambda,
                                   NumericVector integration_constant,
                                   int n_mc,
                                   int n_boot,
                                   bool progress_tf) {

  int indx, i, j;
  NumericMatrix max_dist_uncert(n_boot, n_mc);
  double n;
  double tmp;
  Progress p(n_boot, progress_tf);

  for (j = 0; j < n_boot; j++) {

    for (indx = 0; indx < n_mc; indx++) {

      n = R::rpois(Lambda[j]);
      while (n < 1.0) {
        n = R::rpois(Lambda[j]);
      }

      max_dist_uncert(j, indx) = (sigma[j]/k[j])*(pow(integration_constant[j]*(1 - R::runif(0, 1)), -k[j]) - 1) + mu[j];

      for (i = 1; i < n; i++) {

        tmp = (sigma[j]/k[j])*(pow(integration_constant[j]*(1 - R::runif(0, 1)), -k[j]) - 1) + mu[j];
        if (tmp > max_dist_uncert(j, indx)) {

          max_dist_uncert(j, indx) = tmp;
        }
      }
    }
    p.increment();
  }

  return max_dist_uncert;
}
