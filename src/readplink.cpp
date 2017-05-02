/******************************************************************************/

#include <bigstatsr.h>
#include <fstream>      // std::ofstream

/******************************************************************************/

// [[Rcpp::export]]
bool readbina(const char * filename,
              XPtr<BigMatrix> xpMat,
              const RawMatrix& tab) {
  MatrixAccessor<unsigned char> macc(*xpMat);
  int n = macc.nrow();
  int m = macc.ncol();

  int length = n / 4;
  bool extra = (n > (4 * length));
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

  char c;
  bool is_eof = !(myFile.get(c));
  myFile.close();

  return is_eof;
}

/******************************************************************************/

// [[Rcpp::export]]
void writebina(const char * filename,
               const S4& BM,
               const RawVector& tab,
               const IntegerVector& rowInd,
               const IntegerVector& colInd) {

  XPtr<BigMatrix> xpMat = BM.slot("address");
  RawSubMatAcc macc(*xpMat, rowInd-1, colInd-1, BM.slot("code"));
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
      ind = (macc(i, j) + 4 * macc(i+1, j)) +
        (16 * macc(i+2, j) + 64 * macc(i+3, j));
      buffer[k++] = tab[ind];
    }
    ind = 0; coef = 1;
    for (; i < n; i++) {
      ind += coef * macc(i, j);
      coef *= 4;
    }
    buffer[k] = tab[ind];
    myFile.write(buffer, length); // faster to use (char*)? Nop
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
