/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// Clumping within a distance in bp
// [[Rcpp::export]]
LogicalVector clumping_chr(Environment BM,
                           const IntegerVector& rowInd,
                           const IntegerVector& colInd,
                           const IntegerVector& ordInd,
                           const IntegerVector& pos,
                           const NumericVector& sumX,
                           const NumericVector& denoX,
                           int size,
                           double thr) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);
  int n = macc.nrow();
  int m = macc.ncol();

  double xySum, num, r2;
  int i, j;

  LogicalVector keep(m, false);

  for (int k = 0; k < m; k++) {

    int j0 = ordInd[k] - 1;
    keep[j0] = true;
    int pos_min = pos[j0] - size;
    int pos_max = pos[j0] + size;
    bool not_min = true;
    bool not_max = true;

    for (int l = 1; not_max || not_min; l++) {

      if (not_max) {
        j = j0 + l;
        not_max = (j < m) && (pos[j] <= pos_max);  // within a window..
        if (not_max && keep[j]) {  // look only at already selected ones
          xySum = 0;
          for (i = 0; i < n; i++) {
            xySum += macc(i, j) * macc(i, j0);
          }
          num = xySum - sumX[j] * sumX[j0] / n;
          r2 = num * num / (denoX[j] * denoX[j0]);
          if (r2 > thr) {
            keep[j0] = false;  // prune
            break;
          }
        }
      }

      if (not_min) {
        j = j0 - l;
        not_min = (j >= 0) && (pos[j] >= pos_min);  // within a window..
        if (not_min && keep[j]) {  // look only at already selected ones
          xySum = 0;
          for (i = 0; i < n; i++) {
            xySum += macc(i, j) * macc(i, j0);
          }
          num = xySum - sumX[j] * sumX[j0] / n;
          r2 = num * num / (denoX[j] * denoX[j0]);
          if (r2 > thr) {
            keep[j0] = false;  // prune
            break;
          }
        }
      }
    }
  }

  return keep;
}

/******************************************************************************/
