/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <fstream>
#include "bed-acc.h"

using namespace Rcpp;
using namespace std;

/******************************************************************************/

// [[Rcpp::export]]
void writebina(const char * filename,
               Environment BM,
               const RawVector& tab,
               const IntegerVector& rowInd,
               const IntegerVector& colInd) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  int n = macc.nrow();
  int m = macc.ncol();
  int length = ceil((double)n / 4); // DO NOT USE INTEGERS WITH CEIL

  char *buffer = new char[std::max(3, length) + 1];
  ofstream myFile(filename, ios::out | ios::binary);

  // magic number
  buffer[0] = 108; buffer[1] = 27; buffer[2] = 1;
  myFile.write(buffer, 3);

  int i, j, k, ind, coef;

  for (j = 0; j < m; j++) {
    for (i = 0, k = 0; i <= n-4; i += 4, k++) {
      ind = (macc(i, j) + 4 * macc(i+1, j)) +
        (16 * macc(i+2, j) + 64 * macc(i+3, j));
      buffer[k] = tab[ind];
    }
    ind = 0; coef = 1;
    for (; i < n; i++) {
      ind += coef * macc(i, j);
      coef *= 4;
    }
    buffer[k] = tab[ind];
    myFile.write(buffer, length); // faster to use (char*)? -> no
  }

  myFile.close();
  delete[] buffer;
}

/******************************************************************************/

// for tests
// [[Rcpp::export]]
void testWrite(const RawVector& v, const char * filename) {

  char buffer[256];
  ofstream myFile(filename, ios::out | ios::binary);

  for (int i = 0; i < 256; i++)
    buffer[i] = v[i];

  myFile.write(buffer, 256);
  myFile.close();
}

/******************************************************************************/
