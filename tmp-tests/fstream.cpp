// [[Rcpp::depends(bigmemory, BH)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <fstream>      // std::ifstream
#include <bitset>       // std::bitset

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

