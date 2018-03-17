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

  double w_sum = Rcpp::sum(w), tmp_wsum, tmp_xwsum, tmp_w;

  // beginning
  for (i = 0; i < size; i++) {
    tmp_xwsum = tmp_wsum = 0;
    for (ind_x = i + size, ind_w = w_size - 1; ind_x >= 0; ind_x--, ind_w--) {
      tmp_w = w[ind_w];
      tmp_wsum += tmp_w;
      tmp_xwsum += x[ind_x] * tmp_w;
    }
    res[i] = tmp_xwsum / tmp_wsum;
  }

  // middle
  int lim2 = n - size;
  for (; i < lim2; i++) {
    tmp_xwsum = 0;
    for (ind_x = i - size, ind_w = 0; ind_w < w_size; ind_x++, ind_w++) {
      tmp_xwsum += x[ind_x] * w[ind_w];
    }
    res[i] = tmp_xwsum / w_sum;
  }

  // end
  for (; i < n; i++) {
    tmp_xwsum = tmp_wsum = 0;
    for (ind_x = i - size, ind_w = 0; ind_x < n; ind_x++, ind_w++) {
      tmp_w = w[ind_w];
      tmp_wsum += tmp_w;
      tmp_xwsum += x[ind_x] * tmp_w;
    }
    res[i] = tmp_xwsum / tmp_wsum;
  }

  return res;
}

/******************************************************************************/
