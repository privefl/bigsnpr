/******************************************************************************/

#include <bigstatsr.h>
#include <bigmemory/isna.hpp>
#include <fstream>      // std::ofstream

/******************************************************************************/

// [[Rcpp::export]]
void rawToBigPart(XPtr<BigMatrix> xpMat,
                  const RawVector& source,
                  const RawMatrix& tab,
                  int size, int colOffset,
                  int n, int bsz) {
  MatrixAccessor<unsigned char> macc(*xpMat);

  int i, j, j_off, k, l, c;
  unsigned char t;

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

// [[Rcpp::export]]
void readbina(const char * filename,
              XPtr<BigMatrix> xpMat,
              const RawMatrix& tab) {
  MatrixAccessor<unsigned char> macc(*xpMat);
  int n = macc.nrow();
  int m = macc.ncol();

  int length = n / 4;
  bool extra = (n > 4 * length);
  int lengthExtra = length + extra;

  int i, j, k;
  unsigned char l;

  unsigned char buffer[max(3, lengthExtra)];
  ifstream myFile(filename, ios::in | ios::binary);
  // magic number
  myFile.read((char*)buffer, 3);
  myassert((buffer[0] == 108) && (buffer[1] == 27) && (buffer[2] = 1),
           "Wrong magic number. Aborting..");

  for (j = 0; j < m; j++) {
    // read from bedfile
    myFile.read((char*)buffer, lengthExtra);

    k = 0;
    for (i = 0; i <= n-4; i += 4) {
      l = buffer[k++];
      macc[j][i]   = tab(0, l);
      macc[j][i+1] = tab(1, l);
      macc[j][i+2] = tab(2, l);
      macc[j][i+3] = tab(3, l);
    }
    if (extra) {
      l = buffer[k];
      for (k = 0; i < n; i++, k++) {
        macc[j][i] = tab(k, l);
      }
    }
  }

  myFile.close();
}

/******************************************************************************/

// [[Rcpp::export]]
void writebina(const char * filename,
               XPtr<BigMatrix> xpMat,
               const RawVector& tab) {
  MatrixAccessor<unsigned char> macc(*xpMat);

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
      ind = (macc[j][i] + 4 * macc[j][i+1]) +
        (16 * macc[j][i+2] + 64 * macc[j][i+3]);
      buffer[k++] = tab[ind];
    }
    ind = 0; coef = 1;
    for (; i < n; i++) {
      ind += coef * macc[j][i];
      coef *= 4;
    }
    buffer[k] = tab[ind];
    myFile.write(buffer, length); // faster to use (char*)?
  }

  myFile.close();
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
