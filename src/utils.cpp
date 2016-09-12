// [[Rcpp::depends(RcppEigen, BH, bigmemory)]]
#include <RcppEigen.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;


/******************************************************************************/

// [[Rcpp::export]]
IntegerMatrix mycount(const SEXP pBigMat,
                      const IntegerVector& indCase,
                      const IntegerVector& indControl) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int nCase = indCase.size();
  int nControl = indControl.size();
  int m = xpMat->ncol();

  IntegerMatrix res(6, m);

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < nCase; i++) {
      (res(macc[j][indCase[i]-1], j))++;
    }
    for (int i = 0; i < nControl; i++) {
      (res(macc[j][indControl[i]-1]+3, j))++;
    }
  }

  return(res);
}

/******************************************************************************/

// [[Rcpp::export]]
ListOf<SEXP> mycount2(SEXP pBigMat,
                      const IntegerVector& indCase,
                      const IntegerVector& indControl) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int nCase = indCase.size();
  int nControl = indControl.size();
  int m = xpMat->ncol();

  char tmp;
  int ind;

  IntegerMatrix res(6, m);
  IntegerVector res2(nCase+nControl);

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < nCase; i++) {
      ind = indCase[i]-1;
      tmp = macc[j][ind];
      if (tmp == NA_CHAR) {
        (res2[ind])++;
      } else {
        (res(tmp, j))++;
      }
    }
    for (int i = 0; i < nControl; i++) {
      ind = indControl[i]-1;
      tmp = macc[j][ind];
      if (tmp == NA_CHAR) {
        (res2[ind])++;
      } else {
        (res(tmp+3, j))++;
      }
    }
  }

  return(List::create(_["counts.col"] = res,
                      _["counts.row"] = res2));
}

/******************************************************************************/

// [[Rcpp::export]]
void tcrossprodEigen(SEXP res, const Eigen::Map<Eigen::MatrixXd> bM) {

  XPtr<BigMatrix> xpMatRes(res);

  int n = bM.rows();

  Eigen::Map<Eigen::MatrixXd> bMRes =
    Eigen::Map<Eigen::MatrixXd>((double*)xpMatRes->matrix(), n, n);

  bMRes.selfadjointView<Eigen::Upper>().rankUpdate(bM);

  return;
}

/******************************************************************************/

// [[Rcpp::export]]
void tcrossprodEigen2(SEXP res,
                      const Eigen::Map<Eigen::MatrixXd> X,
                      const Eigen::Map<Eigen::MatrixXd> Y) {

  XPtr<BigMatrix> xpMatRes(res);

  Eigen::Map<Eigen::MatrixXd> bMRes =
    Eigen::Map<Eigen::MatrixXd>((double*)xpMatRes->matrix(),
                                xpMatRes->nrow(),
                                xpMatRes->ncol());

  bMRes += X * Y.transpose();

  return;
}

/******************************************************************************/

// [[Rcpp::export]]
IntegerVector bigcolsumsChar(SEXP pBigMat, const IntegerVector& rowInd) {

  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int n = rowInd.size();
  int m = xpMat->ncol();

  IntegerVector res(m);

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      res[j] += macc[j][rowInd[i]-1];
    }
  }

  return(res);
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector bigcolsumsDouble(SEXP pBigMat) {

  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<double> macc(*xpMat);

  int n = xpMat->nrow();
  int m = xpMat->ncol();

  NumericVector res(m);

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      res[j] += macc[j][i];
    }
  }

  return(res);
}

/******************************************************************************/

// [[Rcpp::export]]
void symCenter(SEXP pBigMat, const NumericVector& means, double mean) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<double> macc(*xpMat);

  int n = xpMat->nrow();

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      macc[j][i] += mean - means[i] - means[j];
    }
  }

  return;
}

/******************************************************************************/

// [[Rcpp::export]]
void colCenter(SEXP pBigMat, const NumericVector& means) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<double> macc(*xpMat);

  int n = xpMat->nrow();
  int m = xpMat->ncol();
  double mj;

  for (int j = 0; j < m; j++) {
    mj = means[j];
    for (int i = 0; i < n; i++) {
      macc[j][i] -= mj;
    }
  }

  return;
}

/******************************************************************************/

// [[Rcpp::export]]
void incrSup(SEXP pBigMat, const NumericMatrix& source) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<double> macc(*xpMat);

  for (int j = 0; j < xpMat->ncol(); j++) {
    for (int i = 0; i <= j; i++) {
      macc[j][i] += source(i,j);
    }
  }

  return;
}

/******************************************************************************/

// [[Rcpp::export]]
void incrAll(SEXP pBigMat, const NumericMatrix& source) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<double> macc(*xpMat);

  for (int j = 0; j < xpMat->ncol(); j++) {
    for (int i = 0; i < xpMat->nrow(); i++) {
      macc[j][i] += source(i,j);
    }
  }

  return;
}

/******************************************************************************/

// [[Rcpp::export]]
void complete(SEXP pBigMat) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<double> macc(*xpMat);

  for (int j = 0; j < xpMat->ncol()-1; j++) {
    for (int i = j+1; i < xpMat->nrow(); i++) {
      macc[j][i] = macc[i][j];
    }
  }

  return;
}

/******************************************************************************/

// [[Rcpp::export]]
NumericMatrix& scaling(NumericMatrix& source,
                       const NumericVector& mean,
                       const NumericVector& sd) {
  int n = source.rows();
  int m = source.cols();

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      source(i,j) -= mean[j];
      source(i,j) /= sd[j];
    }
  }

  return(source);
}

/******************************************************************************/

// [[Rcpp::export]]
NumericMatrix& scaling2(NumericMatrix& source,
                        const NumericVector& intercept,
                        const NumericVector& slope) {
  int n = source.rows();
  int m = source.cols();

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      source(i,j) *= slope[j];
      source(i,j) += intercept[j];
    }
  }

  return(source);
}

/******************************************************************************/

// [[Rcpp::export]]
void rawToBigPart(const IntegerMatrix& source,
                  SEXP pBigMat,
                  int colOffset = 0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int nrows = source.rows();
  int ncols = source.cols();

  for (int j = 0; j < ncols; j++) {
    for (int i = 0; i < nrows; i++) {
      macc[j+colOffset][i] = source(i, j);
    }
  }

  return;
}

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

/******************************************************************************/

template <typename T>
NumericVector R_squared(XPtr<BigMatrix> xpMat,
                        MatrixAccessor<T> macc,
                        const NumericVector& y,
                        const IntegerVector& rowInd,
                        const NumericVector& weights) {
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

// Dispatch function for R_squared
//
// [[Rcpp::export]]
NumericVector R_squared(SEXP pBigMat,
                        const NumericVector& y,
                        const IntegerVector& rowInd,
                        const NumericVector& weights) {
  // First we have to tell Rcpp what class to use for big.matrix objects.
  // This object stores the attributes of the big.matrix object passed to it
  // by R.
  XPtr<BigMatrix> xpMat(pBigMat);

  // To access values in the big.matrix, we need to create a MatrixAccessor
  // object of the appropriate type. Note that in every case we are still
  // returning a NumericVector: this is because big.matrix objects only store
  // numeric values in R, even if their type is set to 'char'. The types
  // simply correspond to the number of bytes used for each element.
  switch(xpMat->matrix_type()) {
  case 1:
    return R_squared(xpMat, MatrixAccessor<char>(*xpMat), y, rowInd, weights);
  case 2:
    return R_squared(xpMat, MatrixAccessor<short>(*xpMat), y, rowInd, weights);
  case 4:
    return R_squared(xpMat, MatrixAccessor<int>(*xpMat), y, rowInd, weights);
  case 8:
    return R_squared(xpMat, MatrixAccessor<double>(*xpMat), y, rowInd, weights);
  default:
    // This case should never be encountered unless the implementation of
    // big.matrix changes, but is necessary to implement shut up compiler
    // warnings.
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

/******************************************************************************/

template <typename T>
NumericMatrix betasRegLin(XPtr<BigMatrix> xpMat,
                          MatrixAccessor<T> macc,
                          const NumericVector& y,
                          const IntegerVector& rowInd,
                          const NumericVector& weights) {
  int n = rowInd.size();
  int m = xpMat->ncol();

  NumericMatrix res(2, m);

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
  double tmp, tmpB;
  double num, denoX;

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
    tmpB = num / denoX;
    res(1, j) = tmpB;
    res(0, j) = (ySum - tmpB * xSum) / wSum;
  }

  return(res);
}

// Dispatch function for betasRegLin
//
// [[Rcpp::export]]
NumericMatrix betasRegLin(SEXP pBigMat,
                          const NumericVector& y,
                          const IntegerVector& rowInd,
                          const NumericVector& weights) {
  // First we have to tell Rcpp what class to use for big.matrix objects.
  // This object stores the attributes of the big.matrix object passed to it
  // by R.
  XPtr<BigMatrix> xpMat(pBigMat);

  // To access values in the big.matrix, we need to create a MatrixAccessor
  // object of the appropriate type. Note that in every case we are still
  // returning a NumericVector: this is because big.matrix objects only store
  // numeric values in R, even if their type is set to 'char'. The types
  // simply correspond to the number of bytes used for each element.
  switch(xpMat->matrix_type()) {
  case 1:
    return betasRegLin(xpMat, MatrixAccessor<char>(*xpMat), y, rowInd, weights);
  case 2:
    return betasRegLin(xpMat, MatrixAccessor<short>(*xpMat), y, rowInd, weights);
  case 4:
    return betasRegLin(xpMat, MatrixAccessor<int>(*xpMat), y, rowInd, weights);
  case 8:
    return betasRegLin(xpMat, MatrixAccessor<double>(*xpMat), y, rowInd, weights);
  default:
    // This case should never be encountered unless the implementation of
    // big.matrix changes, but is necessary to implement shut up compiler
    // warnings.
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

/******************************************************************************/

/* ALREADY IMPLEMENTED IN DEEPCOPY */

// // [[Rcpp::export]]
// void deepcopyPart(SEXP pBigMat,
//                   SEXP pBigMat2,
//                   const IntegerVector& rowInd,
//                   const IntegerVector& colInd) {
//
//   XPtr<BigMatrix> xpMat(pBigMat);
//   MatrixAccessor<char> macc(*xpMat);
//   XPtr<BigMatrix> xpMat2(pBigMat2);
//   MatrixAccessor<char> macc2(*xpMat2);
//
//   int n = rowInd.size();
//   int m = colInd.size();
//
//   int indi, indj;
//
//   for (int j = 0; j < m; j++) {
//     indj = colInd[j] - 1;
//     for (int i = 0; i < n; i++) {
//       indi = rowInd[i] - 1;
//       macc2[j][i] = macc[indj][indi];
//     }
//   }
//
//   return;
// }

/******************************************************************************/
