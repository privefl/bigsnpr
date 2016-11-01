// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;


/******************************************************************************/

// [[Rcpp::export]]
IntegerMatrix rawToBigPart2(SEXP pBigMat,
                            const IntegerVector& rowInd,
                            const IntegerVector& colInd,
                            const IntegerMatrix& tab,
                            const IntegerVector& q,
                            const IntegerVector& r) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<unsigned char> macc(*xpMat);

  int n = rowInd.size();
  int m = colInd.size();

  IntegerVector rows = rowInd - 1;

  IntegerMatrix res(n, m);

  int i, j, k, ind;

  for (j = 0; j < m; j++) {
    k = colInd[j] - 1;
    for (i = 0; i < n; i++) {
      ind = rows[i];
      res(i, j) = tab(r[ind], macc[k][q[ind]]);
    }
  }

  return(res);
}

/******************************************************************************/
