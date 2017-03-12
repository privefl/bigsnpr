// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;

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

// [[Rcpp::export]]
NumericMatrix& toCorr(NumericMatrix& mat,
                      const NumericVector& sqrt_diags) {

  int n = mat.nrow();

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      mat(i, j) /= sqrt_diags[i] * sqrt_diags[j];
    }
  }

  return mat;
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
