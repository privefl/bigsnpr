#include <Rcpp.h>

// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;

// If your compiler is to old, just disable / remove the following line
// [[Rcpp::plugins(cpp11)]]

//' @name R_squared_chr
//' @title Marginally compute R2 with another SNP for each column.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector R_squared_chr(const SEXP pBigMat,
                            const IntegerVector& rowInd,
                            const IntegerVector& colInd,
                            NumericVector colMat0) {

  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int n = rowInd.size();
  int m = colInd.size();
  double nd = (double)n;

  NumericVector res(m);

  double ySum = 0, yySum = 0;
  double tmpY;
  int indi;

  for (int i = 0; i < n; i++) {
    indi = rowInd[i] - 1;
    tmpY = colMat0[indi];
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
      indi = rowInd[i] - 1;
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

// [[Rcpp::export]]
NumericVector R_squared(const SEXP pBigMat,
                        const IntegerVector& y,
                        const IntegerVector& rowInd,
                        const NumericVector& weights) {

  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int n = rowInd.size();
  int m = xpMat->ncol();

  NumericVector res(m);

  double ySum = 0, yySum = 0, wSum = 0;
  double tmpY, tmpW;
  int ind;

  for (int i = 0; i < n; i++) {
    ind = rowInd[i] - 1;
    tmpW = weights[i];
    tmpY = y[ind];
    wSum += tmpW;
    ySum += tmpY * tmpW;
    yySum += tmpY * tmpY * tmpW;
  }
  double denoY = yySum - ySum * ySum / wSum;

  double xSum, xySum, xxSum;
  double tmp;
  double num, denoX;

  for (int j = 0; j < m; j++) {
    xSum = xySum = xxSum = 0;
    for (int i = 0; i < n; i++) {
      ind = rowInd[i] - 1;
      tmpW = weights[i];
      tmp = macc[j][ind];
      xSum += tmp * tmpW;
      xySum += tmp * y[ind] * tmpW;
      xxSum += tmp * tmp * tmpW;
    }
    num = xySum - xSum * ySum / wSum;
    denoX = xxSum - xSum * xSum / wSum;
    res[j] = num * num / (denoX * denoY);
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector betasRegLin(const SEXP pBigMat,
                          const IntegerVector& y,
                          const IntegerVector& rowInd,
                          const NumericVector& weights) {

  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int n = rowInd.size();
  int m = xpMat->ncol();

  NumericVector res(m);

  double ySum = 0, wSum = 0;
  double tmpW;
  int ind;

  for (int i = 0; i < n; i++) {
    ind = rowInd[i] - 1;
    tmpW = weights[i];
    wSum += tmpW;
    ySum += y[ind] * tmpW;
  }

  double xSum, xySum, xxSum;
  double tmp;
  double num, denoX, xSumd;

  for (int j = 0; j < m; j++) {
    xSum = xySum = xxSum = 0;
    for (int i = 0; i < n; i++) {
      ind = rowInd[i] - 1;
      tmp = macc[j][ind];
      tmpW = weights[i];
      xSum += tmp * tmpW;
      xySum += tmp * y[ind] * tmpW;
      xxSum += tmp * tmp * tmpW;
    }
    num = xySum - xSum * ySum / wSum;
    denoX = xxSum - xSum * xSum / wSum;
    res[j] = num / denoX;
  }

  return(res);
}
