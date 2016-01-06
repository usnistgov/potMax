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

// [[Rcpp::export]]
int declusterWithTimeCpp(NumericVector complete_series,
                         NumericVector obs_times,
                         NumericVector y,
                         NumericVector t,
                         double series_mean) {

  int n = complete_series.size();
  int i, indx = 0;
  double curr, prev_high = R_NegInf;
  double curr_time, prev_high_time;

  for (i = 0; i < n; i++) {

    curr = complete_series[i];
    curr_time = obs_times[i];
    if (curr > series_mean) {

      if (curr > prev_high) {

        prev_high = curr;
        prev_high_time = curr_time;
      }
    }
    else {

      if (prev_high > R_NegInf) {

        y[indx] = prev_high;
        t[indx] = prev_high_time;
        indx++;
        prev_high = R_NegInf;
      }
    }
  }

  if (prev_high > R_NegInf) {

    y[indx] = prev_high;
    t[indx] = prev_high_time;
  }

  return(0);
}