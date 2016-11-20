// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector dist1(const IntegerMatrix& mat) {
  int n = mat.ncol();
  //int m = mat.nrow();
  int i, j, dist;
  int c = 0;
  IntegerVector col_i;
  int size = n * (n-1) / 2;
  IntegerVector res(size);
  for (i = 0; i < n; i++) {
    col_i = mat(_, i);
    for (j = 0; j < i; j++) { // only first for now
      dist = sum(abs(col_i - mat(_, j)));
      res[c++] = dist;
    }
  }

  return(res);
}
