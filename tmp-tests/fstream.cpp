// [[Rcpp::depends(bigmemory, BH)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <fstream>      // std::ifstream
#include <bitset>       // std::bitset
#include <bigmemory/isna.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
void readbina(const char * filename,
                       SEXP pBigMat,
                       int length) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int n = xpMat->nrow();
  int m = xpMat->ncol();

  std::ifstream is(filename, std::ifstream::binary);
  std::bitset<8> b;
  char buffer[length];
  int i, j, k, l;

  // magic numbers (already checked in R)
  is.read(buffer, 3);

  // genotypes
  for (j = 0; j < m; j++) {
    // read one SNP at a time
    is.read(buffer,length);
    if (is.fail())
      throw Rcpp::exception("Something went wrong when reading bits");

    // write it in the corresponding column
    i = 0;
    for (k = 0; k < length; k++) {
      b = buffer[k];
      for (l = 0; l < 8 && i < n; l += 2) {
        if (b[l]) {
          macc[j][i++] = b[l+1] ? 2 : NA_CHAR;
        } else {
          macc[j][i++] = b[l+1] ? 1 : 0;
        }
      }
    }
  }

  is.close();
}

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
