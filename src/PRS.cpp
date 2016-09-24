// [[Rcpp::depends(bigmemory, BH)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;


// [[Rcpp::export]]
NumericVector prs2(SEXP pBigMat,
                   const NumericMatrix& odds,
                   const IntegerVector& indTest,
                   const IntegerVector& indCol) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int n = indTest.size();
  int m = indCol.size();

  NumericVector res(n);

  for (int j = 0; j < m; j++) {
    int k = indCol[j]-1;
    for (int i = 0; i < n; i++) {
      res[i] += odds(macc[k][indTest[i]-1], k);
    }
  }

  return(res);
}


// [[Rcpp::export]]
NumericVector prs1(SEXP pBigMat,
                   const NumericVector& betas,
                   const IntegerVector& indTest,
                   const IntegerVector& indCol) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int n = indTest.size();
  int m = indCol.size();

  NumericVector res(n);

  for (int j = 0; j < m; j++) {
    int k = indCol[j]-1;
    for (int i = 0; i < n; i++) {
      res[i] += macc[k][indTest[i]-1] * betas[k];
    }
  }

  return(res);
}
