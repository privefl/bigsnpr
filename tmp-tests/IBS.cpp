/******************************************************************************/

#include "bigsnpr.h"

/******************************************************************************/

// [[Rcpp::export]]
void IBS(XPtr<BigMatrix> xpMat,
         const IntegerVector& code,
         const IntegerMatrix& table,
         XPtr<BigMatrix> xpMat2) {

  // get accessors
  MatrixAccessor<unsigned char> macc(*xpMat);
  int n = macc.nrow();
  int m = macc.ncol();
  MatrixAccessor<int> res(*xpMat2);

  // asserts
  myassert((int)res.nrow()   == 4*n, ERROR_DIM);
  myassert((int)res.ncol()   == n,   ERROR_DIM);
  myassert((int)code.size()  == 256, ERROR_DIM);
  myassert((int)table.nrow() == 4,   ERROR_DIM);
  myassert((int)table.ncol() == 4,   ERROR_DIM);

  int i, j, k;
  unsigned char* macc_j;

  for (j = 0; j < m; j++) {
    macc_j = macc[j];
    for (i = 0; i < n; i++) {
      for (k = 0; k < i; k++) {
        (res[i][4 * k + table(code[macc_j[i]],
                              code[macc_j[k]])])++;
      }
    }
  }
}

/******************************************************************************/


// [[Rcpp::export]]
arma::Cube<int> IBS2(XPtr<BigMatrix> xpMat,
                const IntegerVector& code,
                const IntegerMatrix& table,
                IntegerVector& res) {

  // get accessors
  MatrixAccessor<unsigned char> macc(*xpMat);
  int n = macc.nrow();
  int m = macc.ncol();

  IntegerVector dimCube = res.attr("dim");
  arma::Cube<int> res_(res.begin(), dimCube[0], dimCube[1], dimCube[2]);

  // asserts
  myassert((int)res.size()   == 4*n*n, ERROR_DIM);
  myassert((int)code.size()  == 256,   ERROR_DIM);
  myassert((int)table.nrow() == 4,     ERROR_DIM);
  myassert((int)table.ncol() == 4,     ERROR_DIM);

  int i, j, k, ind;
  unsigned char* macc_j;

  for (j = 0; j < m; j++) {
    macc_j = macc[j];
    for (i = 0; i < n; i++) {
      for (k = 0; k < i; k++) {
        ind = table(code[macc_j[k]], code[macc_j[i]]);
        (res_(ind, k, i))++;
      }
    }
  }

  return res_;
}

/******************************************************************************/
