#include <Rcpp.h>
using namespace Rcpp;

#include <RProgress.h>

// [[Rcpp::export]]
NumericVector gumbelMaxDistCpp(double mu, double sigma,
                               double Lambda,
                               double integration_constant,
                               int n_mc, bool progress_tf) {

  NumericVector max_dist(n_mc);
  int indx, i;
  double n;
  RProgress::RProgress p("|:bar| :percent ~ :eta",
                         n_mc,
                         0.6*Rf_GetOptionWidth(),
                         '+',
                         ' ',
                         false,
                         0.2);
  double tmp;

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

    if (progress_tf) {
      p.tick();
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
                                     int n_boot,
                                     bool progress_tf) {

  int indx, i, j;
  NumericMatrix max_dist_uncert(n_boot, n_mc);
  double n;
  double tmp, tmp_max;
  RProgress::RProgress p("|:bar| :percent ~ :eta",
                         n_boot,
                         0.6*Rf_GetOptionWidth(),
                         '+',
                         ' ',
                         false,
                         0.2);

  for (j = 0; j < n_boot; j++) {

    for (indx = 0; indx < n_mc; indx++) {

      n = R::rpois(Lambda[j]);
      while (n < 1.0) {
        n = R::rpois(Lambda[j]);
      }

      tmp_max = mu[j] - sigma[j]*(log(integration_constant[j]) + log(1 - R::runif(0, 1)));

      for (i = 1; i < n; i++) {

        tmp = mu[j] - sigma[j]*(log(integration_constant[j]) + log(1 - R::runif(0, 1)));
        if (tmp > tmp_max) {

          tmp_max = tmp;
        }
      }
      max_dist_uncert(j, indx) = tmp_max;
    }
    if(progress_tf) {
      p.tick();
    }
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
  RProgress::RProgress p("|:bar| :percent ~ :eta",
                         n_mc,
                         0.6*Rf_GetOptionWidth(),
                         '+',
                         ' ',
                         false,
                         0.2);

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
    if (progress_tf) {
      p.tick();
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
                                   int n_boot,
                                   bool progress_tf) {

  int indx, i, j;
  NumericMatrix max_dist_uncert(n_boot, n_mc);
  double n;
  double tmp, tmp_max;
  RProgress::RProgress p("|:bar| :percent ~ :eta",
                         n_boot,
                         0.6*Rf_GetOptionWidth(),
                         '+',
                         ' ',
                         false,
                         0.2);

  for (j = 0; j < n_boot; j++) {

    for (indx = 0; indx < n_mc; indx++) {

      n = R::rpois(Lambda[j]);
      while (n < 1.0) {
        n = R::rpois(Lambda[j]);
      }

      tmp_max = (sigma[j]/k[j])*(pow(integration_constant[j]*(1 - R::runif(0, 1)), -k[j]) - 1) + mu[j];

      for (i = 1; i < n; i++) {

        tmp = (sigma[j]/k[j])*(pow(integration_constant[j]*(1 - R::runif(0, 1)), -k[j]) - 1) + mu[j];
        if (tmp > tmp_max) {

          tmp_max = tmp;
        }
      }
      max_dist_uncert(j, indx) = tmp_max;
    }
    if (progress_tf) {
      p.tick();
    }
  }

  return max_dist_uncert;
}
