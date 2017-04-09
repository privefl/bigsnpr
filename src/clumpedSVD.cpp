/******************************************************************************/

#include "bigsnpr.h"

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

// Clumping in a restricted number of SNPs
// [[Rcpp::export]]
LogicalVector local_clumping(const S4& BM,
                             const IntegerVector& rowInd,
                             const IntegerVector& colInd,
                             const IntegerVector& ordInd,
                             const NumericVector& sumX,
                             const NumericVector& denoX,
                             double thr) {

  XPtr<BigMatrix> xpMat = BM.slot("address");
  RawSubMatAcc macc(*xpMat, rowInd-1, colInd-1, BM.slot("code"));

  int n = macc.nrow();
  int m = macc.ncol();

  double xySum, num, r2;
  int i, j, j0, k;

  LogicalVector remain(m, true); // init with all true
  LogicalVector keep(m); // init with all false

  for (k = 0; k < m; k++) {
    j0 = ordInd[k] - 1;
    if (remain[j0]) { // if already excluded, goto next
      remain[j0] = false;
      keep[j0] = true;
      for (j = 0; j < m; j++) {
        if (remain[j]) { // if already excluded, goto next
          xySum = 0;
          for (i = 0; i < n; i++) {
            xySum += macc(i, j) * macc(i, j0);
          }
          num = xySum - sumX[j] * sumX[j0] / n;
          r2 = num * num / (denoX[j] * denoX[j0]);
          if (r2 > thr) remain[j] = false; // prune
        }
      }
    }
  }

  return keep;
}

/******************************************************************************/
