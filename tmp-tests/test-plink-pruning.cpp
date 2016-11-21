// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;


/******************************************************************************/

// [[Rcpp::export]]
LogicalVector& R_squared_chr2(SEXP pBigMat,
                              LogicalVector& keep,
                              const NumericVector& mafX,
                              const NumericVector& sumX,
                              const NumericVector& sdX,
                              int size,
                              double thr) {
  // Assert that keep[j] == TRUE
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int n = xpMat->nrow();
  double nd = (double)n;
  int m = xpMat->ncol();
  double xySum, num, r2;

  int j0, j, i;

  for (j0 = 1; j0 < size; j0++) {
    if (keep[j0]) { // if already pruned, goto next
      for (j = 0; j < j0; j++) {
        if (keep[j]) { // if already pruned, goto next
          xySum = 0;
          for (i = 0; i < n; i++) {
            xySum += macc[j][i] * macc[j0][i];
          }
          num = xySum - sumX[j] * sumX[j0] / nd;
          r2 = num * num / (sdX[j] * sdX[j0]);
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

  for(j0 = size; j0 < m; j0++) {
    if (keep[j0]) { // if already pruned, goto next
      for (j = j0 - size + 1; j < j0; j++) {
        if (keep[j]) { // if already pruned, goto next
          xySum = 0;
          for (i = 0; i < n; i++) {
            xySum += macc[j][i] * macc[j0][i];
          }
          num = xySum - sumX[j] * sumX[j0] / nd;
          r2 = num * num / (sdX[j] * sdX[j0]);
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

  return(keep);
}

/******************************************************************************/
