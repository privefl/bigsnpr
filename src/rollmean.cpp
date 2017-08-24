/******************************************************************************/

#include <Rcpp.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
NumericVector roll_mean(const NumericVector& x,
                        const NumericVector& w) {

  int n = x.size();
  int w_size = w.size();
  int size = (w_size - 1) / 2;

  NumericVector res(n);
  int i, ind_x, ind_w;

  double w_sum = Rcpp::sum(w);
  double tmp_wsum;

  for (i = 0; i < n; i++) {
    if ((i - size) < 0) { // beginning
      tmp_wsum = 0;
      for (ind_x = i + size, ind_w = w_size - 1; ind_x >= 0; ind_x--, ind_w--) {
        res[i] += x[ind_x] * w[ind_w];
        tmp_wsum += w[ind_w];
      }
      res[i] /= tmp_wsum;
    } else if ((i + size) >= n) { // end
      tmp_wsum = 0;
      for (ind_x = i - size, ind_w = 0; ind_x < n; ind_x++, ind_w++) {
        res[i] += x[ind_x] * w[ind_w];
        tmp_wsum += w[ind_w];
      }
      res[i] /= tmp_wsum;
    } else { // middle
      for (ind_x = i - size, ind_w = 0; ind_w < w_size; ind_x++, ind_w++) {
        res[i] += x[ind_x] * w[ind_w];
      }
      res[i] /= w_sum;
    }
  }

  return res;
}

/******************************************************************************/
