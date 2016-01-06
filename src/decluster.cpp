#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int declusterCpp(NumericVector complete_series,
                 NumericVector y,
                 double series_mean) {

  int n = complete_series.size();
  int i, indx = 0;
  double curr, prev_high = R_NegInf;

  for (i = 0; i < n; i++) {

    curr = complete_series[i];
    if (curr > series_mean) {

      if (curr > prev_high) {

        prev_high = curr;
      }
    }
    else {

      if (prev_high > R_NegInf) {

        y[indx] = prev_high;
        indx++;
        prev_high = R_NegInf;
      }
    }
  }

  if (prev_high > R_NegInf) {

    y[indx] = prev_high;
  }

  return(0);
}
