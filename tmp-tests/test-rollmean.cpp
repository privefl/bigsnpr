/******************************************************************************/

#include <Rcpp.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
NumericVector roll_mean0(const NumericVector& x,
                         const NumericVector& w) {

  int n = x.size();
  int w_size = w.size();
  int size = (w_size - 1) / 2;

  NumericVector res(n);
  int i, ind_x, ind_w;

  double w_sum = Rcpp::sum(w);
  double tmp_wsum;

  for (i = 0; i < n; i++) {
    // Rcout << "ONE" << std::endl;
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

// [[Rcpp::export]]
NumericVector roll_mean(const NumericVector& x,
                        const NumericVector& w) {

  int n = x.size();
  int w_size = w.size();
  int size = (w_size - 1) / 2;

  NumericVector res(n);
  int i, ind_x, ind_w;

  double w_sum = Rcpp::sum(w), tmp_wsum;

  // beginning
  for (i = 0; i < size; i++) {
    tmp_wsum = 0;
    for (ind_x = i + size, ind_w = w_size - 1; ind_x >= 0; ind_x--, ind_w--) {
      res[i] += x[ind_x] * w[ind_w];
      tmp_wsum += w[ind_w];
    }
    res[i] /= tmp_wsum;
  }

  // middle
  int lim2 = n - size;
  for (; i < lim2; i++) {
    for (ind_x = i - size, ind_w = 0; ind_w < w_size; ind_x++, ind_w++) {
      res[i] += x[ind_x] * w[ind_w];
    }
    res[i] /= w_sum;
  }

  // end
  for (; i < n; i++) {
    tmp_wsum = 0;
    for (ind_x = i - size, ind_w = 0; ind_x < n; ind_x++, ind_w++) {
      res[i] += x[ind_x] * w[ind_w];
      tmp_wsum += w[ind_w];
    }
    res[i] /= tmp_wsum;
  }

  return res;
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector roll_mean2(const NumericVector& x,
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



/*** R
for (add_x in 0:1) {
  for (add_w in 0:1) {
    x <- runif(1e5 + add_x)
    w <- runif(50 + add_w)
    stopifnot(isTRUE(all.equal(roll_mean0(x, w), roll_mean(x, w))))
    stopifnot(isTRUE(all.equal(roll_mean0(x, w), roll_mean2(x, w))))
    print(microbenchmark::microbenchmark(
      roll_mean0(x, w),
      roll_mean(x, w),
      roll_mean2(x, w)
    ))
  }
}
*/
