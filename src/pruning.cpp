// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;


/******************************************************************************/

// [[Rcpp::export]]
NumericVector R_squared_chr(SEXP pBigMat,
                            const IntegerVector& rowInd,
                            const IntegerVector& colInd,
                            const NumericVector& colMat0) {

  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int n = rowInd.size();
  int m = colInd.size();
  double nd = (double)n;

  // indices begin at 1 in R and 0 in C++
  IntegerVector trains = rowInd - 1;

  NumericVector res(m);

  double ySum = 0, yySum = 0;
  double tmpY;
  int indi;

  for (int i = 0; i < n; i++) {
    tmpY = colMat0[trains[i]];
    ySum += tmpY;
    yySum += tmpY * tmpY;
  }
  double denoY = yySum - ySum * ySum / nd;

  int xSum, xySum, xxSum;
  char tmp;
  int indj;
  double num, denoX, xSumd;

  for (int j = 0; j < m; j++) {
    indj = colInd[j] - 1;
    xSum = xySum = xxSum = 0;
    for (int i = 0; i < n; i++) {
      indi = trains[i];
      tmp = macc[indj][indi];
      xSum += tmp;
      xySum += tmp * colMat0[indi];
      xxSum += tmp * tmp;
    }
    xSumd = (double)xSum;
    num = (double)xySum - xSumd * ySum / nd;
    denoX = (double)xxSum - xSumd * xSumd / nd;
    res[j] = num * num / (denoX * denoY);
  }

  return(res);
}

/******************************************************************************/

// [[Rcpp::export]]
LogicalVector& R_squared_chr2(SEXP pBigMat,
                              const IntegerVector& rowInd,
                              LogicalVector& keep,
                              const NumericVector& mafX,
                              const NumericVector& sumX,
                              const NumericVector& denoX,
                              int size,
                              double thr) {
  // Assert that keep[j] == TRUE
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  // indices begin at 1 in R and 0 in C++
  IntegerVector trains = rowInd - 1;

  int n = rowInd.size();
  double nd = (double)n;
  int m = xpMat->ncol();
  double xySum, num, r2;

  int j0, j, i, ind_i;

  for (j0 = 1; j0 < size; j0++) {
    if (keep[j0]) { // if already pruned, goto next
      for (j = 0; j < j0; j++) {
        if (keep[j]) { // if already pruned, goto next
          xySum = 0;
          for (i = 0; i < n; i++) {
            ind_i = trains[i];
            xySum += macc[j][ind_i] * macc[j0][ind_i];
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

  for(j0 = size; j0 < m; j0++) {
    if (keep[j0]) { // if already pruned, goto next
      for (j = j0 - size + 1; j < j0; j++) {
        if (keep[j]) { // if already pruned, goto next
          xySum = 0;
          for (i = 0; i < n; i++) {
            ind_i = trains[i];
            xySum += macc[j][ind_i] * macc[j0][ind_i];
          }
          num = xySum - sumX[j] * sumX[j0] / nd;
          r2 = num * num / (denoX[j] * denoX[j0]);
          if (r2 > thr) { // prune one of them
            if (mafX[j] < mafX[j0]) { // prune the one with smaller maf
              keep[j] = false;
            } else {
              keep[j0] = false;
              break;
            }
          }
        }
      }
    }
  }

  return(keep);
}

/******************************************************************************/
