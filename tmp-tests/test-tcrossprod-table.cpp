// [[Rcpp::depends(bigmemory, BH, RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix& tabcrossprod(const IntegerMatrix& mat,
                            int n, int m,
                            NumericMatrix& res,
                            const arma::cube& tab) {
  double tmp;
  int i, j, k;

  for (j = 0; j < n; j++) {
    for (i = j; i < n; i++) { // lower tri
      tmp = 0;
      for (k = 0; k < m; k++) {
        tmp += tab(mat(i, k), mat(j, k), k);
      }
      res(i, j) += tmp;
    }
  }

  return(res);
}


