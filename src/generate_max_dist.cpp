#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

// [[Rcpp::export]]
NumericVector gumbelMaxDistCpp(double mu, double sigma,
                               double Lambda,
                               double integration_constant,
                               int n_mc) {

  NumericVector max_dist(n_mc);
  int indx = 0, i;
  double n;
  double tmp;
  Progress p(n_mc, true);

  while (indx < n_mc) {

    n = R::rpois(Lambda);

    if (n > 0) {

      max_dist[indx] = mu - sigma*(log(integration_constant) + log(1 - R::runif(0, 1)));

      for (i = 1; i < n; i++) {

        tmp = mu - sigma*(log(integration_constant) + log(1 - R::runif(0, 1)));
        if (tmp > max_dist[indx]) {

          max_dist[indx] = tmp;
        }
      }

      indx++;
      p.increment();
    }
  }

  return max_dist;
}

// [[Rcpp::export]]
NumericMatrix gumbelMaxDistUncertCpp(NumericVector mu,
                                     NumericVector sigma,
                                     NumericVector Lambda,
                                     NumericVector integration_constant,
                                     int n_mc,
                                     int n_boot) {

  int indx, i, j;
  NumericMatrix max_dist_uncert(n_boot, n_mc);
  double n;
  double tmp;
  Progress p(n_boot, true);

  for (j = 0; j < n_boot; j++) {

    indx = 0;

    while (indx < n_mc) {

      n = R::rpois(Lambda[j]);

      if (n > 0) {

        max_dist_uncert(j, indx) = mu[j] - sigma[j]*(log(integration_constant[j]) + log(1 - R::runif(0, 1)));

        for (i = 1; i < n; i++) {

          tmp = mu[j] - sigma[j]*(log(integration_constant[j]) + log(1 - R::runif(0, 1)));
          if (tmp > max_dist_uncert(j, indx)) {

            max_dist_uncert(j, indx) = tmp;
          }
        }

        indx++;
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
                             int n_mc) {

  NumericVector max_dist(n_mc);
  int indx = 0, i;
  double n;
  double tmp;
  Progress p(n_mc, true);

  while (indx < n_mc) {

    n = R::rpois(Lambda);

    if (n > 0) {

      max_dist[indx] = (sigma/k)*(pow(integration_constant*(1 - R::runif(0, 1)), -k) - 1) + mu;

      for (i = 1; i < n; i++) {

        tmp = (sigma/k)*(pow(integration_constant*(1 - R::runif(0, 1)), -k) - 1) + mu;
        if (tmp > max_dist[indx]) {

          max_dist[indx] = tmp;
        }
      }

      indx++;
      p.increment();
    }
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
                                   int n_boot) {

  int indx, i, j;
  NumericMatrix max_dist_uncert(n_boot, n_mc);
  double n;
  double tmp;
  Progress p(n_boot, true);

  for (j = 0; j < n_boot; j++) {

    indx = 0;

    while (indx < n_mc) {

      n = R::rpois(Lambda[j]);

      if (n > 0) {

        max_dist_uncert(j, indx) = (sigma[j]/k[j])*(pow(integration_constant[j]*(1 - R::runif(0, 1)), -k[j]) - 1) + mu[j];

        for (i = 1; i < n; i++) {

          tmp = (sigma[j]/k[j])*(pow(integration_constant[j]*(1 - R::runif(0, 1)), -k[j]) - 1) + mu[j];
          if (tmp > max_dist_uncert(j, indx)) {

            max_dist_uncert(j, indx) = tmp;
          }
        }

        indx++;
      }
    }
    p.increment();
  }

  return max_dist_uncert;
}
