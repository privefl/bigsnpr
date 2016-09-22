// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;


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
