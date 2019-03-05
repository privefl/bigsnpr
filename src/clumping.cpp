/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

using namespace Rcpp;

/******************************************************************************/

// Clumping within a distance in bp
// [[Rcpp::export]]
LogicalVector clumping2(Environment BM,
                        const IntegerVector& rowInd,
                        const IntegerVector& colInd,
                        const IntegerVector& ordInd,
                        const IntegerVector& pos,
                        const NumericVector& sumX,
                        const NumericVector& denoX,
                        int size,
                        double thr) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd - 1, colInd - 1, BM["code256"]);
  int n = macc.nrow();
  int m = macc.ncol();

  double xySum, num, r2;
  int i, j, j0, k, l, pos_min, pos_max;
  bool cont1, cont2;

  LogicalVector keep(m, false);

  for (k = 0; k < m; k++) {

    j0 = ordInd[k] - 1;
    keep[j0] = true;
    cont1 = cont2 = true;
    pos_min = pos[j0] - size;
    pos_max = pos[j0] + size;

    for (l = 1; cont1 || cont2; l++) {

      if (cont1) {
        j = j0 + l;
        cont1 = (j < m) && (pos[j] <= pos_max);  // within a window..
        if (cont1 && keep[j]) {  // look only at already selected ones
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

      if (cont2) {
        j = j0 - l;
        cont2 = (j >= 0) && (pos[j] >= pos_min);  // within a window..
        if (cont2 && keep[j]) {  // look only at already selected ones
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
