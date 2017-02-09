// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;


/******************************************************************************/

// [[Rcpp::export]]
void doubleBM(XPtr<BigMatrix> xpMat, XPtr<BigMatrix> xpMat2) {
  MatrixAccessor<char> macc(*xpMat);
  MatrixAccessor<char> macc2(*xpMat2);

  int n = macc.nrow();
  int m = macc.ncol();

  int i, j, j2;
  char tmp;

  for (j = j2 = 0; j < m; j++, j2 += 2) {
    for (i = 0; i < n; i++) {
      tmp = macc[j][i];
      if (tmp == 0) {
        macc2[j2][i] = macc2[j2+1][i] = 0;
      } else if (tmp == 1) {
        macc2[j2][i] = 0;
        macc2[j2+1][i] = 1;
      } else if (tmp == 2) {
        macc2[j2][i] = macc2[j2+1][i] = 1;
      } else {
        throw Rcpp::exception("Your big.matrix should have only Os, 1s or 2s");
      }
    }
  }
}

/******************************************************************************/
