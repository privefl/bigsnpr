// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/isna.hpp>
#include <fstream>      // std::ofstream

using namespace Rcpp;


/******************************************************************************/

// [[Rcpp::export]]
void rawToBigPart(XPtr<BigMatrix> xpMat,
                  const RawVector& source,
                  const IntegerMatrix& tab,
                  int size, int colOffset,
                  int n, int bsz) {
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

inline char replaceNA(char x) {
  return(isna(x) ? 3 : x);
}

// [[Rcpp::export]]
void writebina(const char * filename,
               XPtr<BigMatrix> xpMat,
               const RawVector& tab) {
  MatrixAccessor<char> macc(*xpMat);

  int n = macc.nrow();
  int m = macc.ncol();
  int length = ceil((double)n / 4); // DO NOT USE INTEGERS WITH CEIL

  char buffer[max(3, length)];
  ofstream myFile(filename, ios::out | ios::binary);

  // magic number
  buffer[0] = 108; buffer[1] = 27; buffer[2] = 1;
  myFile.write(buffer, 3);

  int i, j, k, ind, coef;

  for (j = 0; j < m; j++) {
    k = 0;
    for (i = 0; i <= n-4; i += 4) {
      ind = (replaceNA(macc[j][i]) + 4 * replaceNA(macc[j][i+1])) +
        (16 * replaceNA(macc[j][i+2]) + 64 * replaceNA(macc[j][i+3]));
      buffer[k++] = tab[ind];
    }
    ind = 0; coef = 1;
    for (; i < n; i++) {
      ind += coef * replaceNA(macc[j][i]);
      coef *= 4;
    }
    buffer[k] = tab[ind];
    myFile.write(buffer, length);
  }

  myFile.close();
}

/******************************************************************************/
