// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;


/******************************************************************************/

// [[Rcpp::export]]
Rcpp::NumericVector R_squared_chr(SEXP pBigMat,
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
