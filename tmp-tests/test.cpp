// [[Rcpp::depends(RcppEigen, BH, bigmemory)]]
#include <RcppEigen.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
void mycount(const SEXP pBigMat,
                      const IntegerVector& indCase,
                      const IntegerVector& indControl) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int n = xpMat->nrow();
  int m = xpMat->ncol();

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      if (Rcpp::is_na(macc[j][i]))
        printf("NA ");
    }
  }

  return;
}
