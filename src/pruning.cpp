/******************************************************************************/

#include "bigsnpr.h"

/******************************************************************************/

// [[Rcpp::export]]
LogicalVector clumping(XPtr<BigMatrix> xpMat,
                       const IntegerVector& rowInd,
                       const IntegerVector& colInd,
                       const IntegerVector& ordInd,
                       LogicalVector& remain,
                       const NumericVector& sumX,
                       const NumericVector& denoX,
                       int size,
                       double thr) {
  SubMatrixAccessor<char> macc(*xpMat, rowInd-1, colInd-1);

  int n = macc.nrow();
  int m = macc.ncol();
  double nd = (double)n;

  double xySum, num, r2;
  int i, j, j0, k;

  LogicalVector keep(m); // init with all false

  for (k = 0; k < m; k++) {
    j0 = ordInd[k] - 1;
    if (remain[j0]) { // if already excluded, goto next
      for (j = max(0, j0 - size); j < min(m, j0 + size); j++) {
        if (remain[j]) { // if already excluded, goto next
          xySum = 0;
          for (i = 0; i < n; i++) {
            xySum += macc(i, j) * macc(i, j0);
          }
          num = xySum - sumX[j] * sumX[j0] / nd;
          r2 = num * num / (denoX[j] * denoX[j0]);
          if (r2 > thr) remain[j] = false; // prune
        }
      }
      keep[j0] = true;
      remain[j0] = false;
    }
  }

  return keep;
}

/******************************************************************************/

// Pruning within a distance in number of SNPs
// [[Rcpp::export]]
LogicalVector& pruning(XPtr<BigMatrix> xpMat,
                       const IntegerVector& rowInd,
                       const IntegerVector& colInd,
                       LogicalVector& keep,
                       const NumericVector& mafX,
                       const NumericVector& sumX,
                       const NumericVector& denoX,
                       int size,
                       double thr) {
  // Assert that keep[j] == TRUE
  SubMatrixAccessor<char> macc(*xpMat, rowInd-1, colInd-1);

  int n = macc.nrow();
  double nd = (double)n;
  int m = macc.ncol();
  double xySum, num, r2;

  int j0, j, i, j_max;

  for (j0 = 0; j0 < m; j0++) {
    if (keep[j0]) { // if already excluded, goto next
      j_max = min(j0 + size, m);
      for (j = j0 + 1; j < j_max; j++) {
        if (keep[j]) { // if already excluded, goto next
          xySum = 0;
          for (i = 0; i < n; i++) {
            xySum += macc(i, j) * macc(i, j0);
          }
          num = xySum - sumX[j] * sumX[j0] / nd;
          r2 = num * num / (denoX[j] * denoX[j0]);
          if (r2 > thr) { // prune one of them
            if (mafX[j0] < mafX[j]) { // prune the one with smaller maf
              keep[j0] = false;
              break;
            } else {
              keep[j] = false;
            }
          }
        }
      }
    }
  }

  return keep;
}

// Pruning within a distance in kb
// [[Rcpp::export]]
LogicalVector& pruning2(XPtr<BigMatrix> xpMat,
                        const IntegerVector& rowInd,
                        const IntegerVector& colInd,
                        LogicalVector& keep,
                        const IntegerVector& pos,
                        const NumericVector& mafX,
                        const NumericVector& sumX,
                        const NumericVector& denoX,
                        int size,
                        double thr) {
  // Assert that keep[j] == TRUE
  SubMatrixAccessor<char> macc(*xpMat, rowInd-1, colInd-1);

  int n = macc.nrow();
  double nd = (double)n;
  int m = macc.ncol();
  double xySum, num, r2;
  int pos_max;

  int j0, j, i;

  for (j0 = 0; j0 < m; j0++) {
    if (keep[j0]) { // if already excluded, goto next
      pos_max = pos[j0] + size;
      for (j = j0 + 1; pos[j] <= pos_max; j++) { // pos[m] == 0 -> break
        if (keep[j]) { // if already excluded, goto next
          xySum = 0;
          for (i = 0; i < n; i++) {
            xySum += macc(i, j) * macc(i, j0);
          }
          num = xySum - sumX[j] * sumX[j0] / nd;
          r2 = num * num / (denoX[j] * denoX[j0]);
          if (r2 > thr) { // prune one of them
            if (mafX[j0] < mafX[j]) { // prune the one with smaller maf
              keep[j0] = false;
              break;
            } else {
              keep[j] = false;
            }
          }
        }
      }
    }
  }

  return keep;
}

/******************************************************************************/
