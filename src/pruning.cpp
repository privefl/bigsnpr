/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

using namespace Rcpp;

/******************************************************************************/

// Clumping within a distance in number of SNPs
// [[Rcpp::export]]
LogicalVector clumping(Environment BM,
                       const IntegerVector& rowInd,
                       const IntegerVector& colInd,
                       const IntegerVector& ordInd,
                       LogicalVector& remain,
                       const NumericVector& sumX,
                       const NumericVector& denoX,
                       int size,
                       double thr) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd - 1, colInd - 1, BM["code256"]);
  int n = macc.nrow();
  int m = macc.ncol();

  double xySum, num, r2;
  int i, j, j0, k, j_min, j_max;

  LogicalVector keep(m); // init with all false

  for (k = 0; k < m; k++) {
    j0 = ordInd[k] - 1;
    if (remain[j0]) { // if already excluded, goto next
      remain[j0] = false;
      keep[j0] = true;
      j_min = std::max(0, j0 - size);
      j_max = std::min(m, j0 + size + 1);
      for (j = j_min; j < j_max; j++) {
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

// Clumping within a distance in bp
// [[Rcpp::export]]
LogicalVector clumping2(Environment BM,
                        const IntegerVector& rowInd,
                        const IntegerVector& colInd,
                        const IntegerVector& ordInd,
                        LogicalVector& remain,
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
  int i, j, j0, k, pos_min, pos_max;

  LogicalVector keep(m); // init with all false

  for (k = 0; k < m; k++) {
    j0 = ordInd[k] - 1;
    if (remain[j0]) {
      remain[j0] = false;
      keep[j0] = true;
      pos_min = pos[j0] - size;
      pos_max = pos[j0] + size;
      for (j = 0; pos[j] <= pos_max; j++) { // pos[m] == MAX -> break
        if (remain[j] && (pos[j] >= pos_min)) {
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

// Pruning within a distance in number of SNPs
// [[Rcpp::export]]
LogicalVector& pruning(Environment BM,
                       const IntegerVector& rowInd,
                       const IntegerVector& colInd,
                       LogicalVector& keep,
                       const NumericVector& mafX,
                       const NumericVector& sumX,
                       const NumericVector& denoX,
                       int size,
                       double thr) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd - 1, colInd - 1, BM["code256"]);
  int n = macc.nrow();
  int m = macc.ncol();

  double xySum, num, r2;
  int j0, j, i, j_max;

  for (j0 = 0; j0 < m; j0++) {
    if (keep[j0]) { // if already excluded, goto next
      j_max = std::min(j0 + size + 1, m);
      for (j = j0 + 1; j < j_max; j++) {
        if (keep[j]) { // if already excluded, goto next
          xySum = 0;
          for (i = 0; i < n; i++) {
            xySum += macc(i, j) * macc(i, j0);
          }
          num = xySum - sumX[j] * sumX[j0] / n;
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

// Pruning within a distance in bp
// [[Rcpp::export]]
LogicalVector& pruning2(Environment BM,
                        const IntegerVector& rowInd,
                        const IntegerVector& colInd,
                        LogicalVector& keep,
                        const IntegerVector& pos,
                        const NumericVector& mafX,
                        const NumericVector& sumX,
                        const NumericVector& denoX,
                        int size,
                        double thr) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd - 1, colInd - 1, BM["code256"]);
  int n = macc.nrow();
  int m = macc.ncol();

  double xySum, num, r2;
  int j0, j, i, pos_max;

  for (j0 = 0; j0 < m; j0++) {
    if (keep[j0]) {
      pos_max = pos[j0] + size;
      for (j = j0 + 1; pos[j] <= pos_max; j++) { // pos[m] == MAX -> break
        if (keep[j]) {
          xySum = 0;
          for (i = 0; i < n; i++) {
            xySum += macc(i, j) * macc(i, j0);
          }
          num = xySum - sumX[j] * sumX[j0] / n;
          r2 = num * num / (denoX[j] * denoX[j0]);
          if (r2 > thr) {
            if (mafX[j0] < mafX[j]) {
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
