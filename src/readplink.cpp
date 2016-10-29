// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;


/******************************************************************************/

// [[Rcpp::export]]
void rawToBigPart(SEXP pBigMat,
                  const RawVector& source,
                  const IntegerMatrix& tab,
                  int size, int colOffset,
                  int n, int bsz) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int i, j, j_off, k, l, t, c = 0;

  for (j = 0; j < size; j++) {
    j_off = j + colOffset;
    i = 0;
    for (k = 0; k < bsz; k++) {
      t = source[c++];
      for (l = 0; l < 4 && i < n; l++) {
        macc[j_off][i++] = tab(l, t);
      }
    }
  }
}

/******************************************************************************/
