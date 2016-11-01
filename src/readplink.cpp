// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/isna.hpp>
#include <fstream>      // std::ofstream

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

  int i, j, j_off, k, l, t, c;

  c = 0;
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

char replaceNA(char x) {
  return(isna(x) ? 3 : x);
}

// [[Rcpp::export]]
void writebina(const char * filename,
               SEXP pBigMat,
               const RawVector& tab) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int n = xpMat->nrow();
  int q = n / 4;
  int r = n % 4;
  int m = xpMat->ncol();

  int length = q + (r > 0);
  char buffer[length];
  ofstream myFile(filename, ios::out | ios::binary);
  buffer[0] = 108; buffer[1] = 27; buffer[2] = 1;
  myFile.write(buffer, 3);

  int i, j, k, ind;

  for (j = 0; j < m; j++) {
    i = 0;
    for (k = 0; k < q; k++) {
      ind = replaceNA(macc[j][i++]);
      ind += 4 * replaceNA(macc[j][i++]);
      ind += 16 * replaceNA(macc[j][i++]);
      ind += 64 * replaceNA(macc[j][i++]);
      buffer[k] = tab[ind];
    }
    if (r > 0) {
      ind = replaceNA(macc[j][i++]);
      if (r > 1) {
        ind += 4 * replaceNA(macc[j][i++]);
        if (r > 2) {
          ind += 16 * replaceNA(macc[j][i++]);
        }
      }
      buffer[q] = tab[ind];
    }
    myFile.write(buffer, length);
  }
  myFile.close();
}

/******************************************************************************/
