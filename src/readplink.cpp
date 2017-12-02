/******************************************************************************/

#include <RcppArmadillo.h>
#include <bigstatsr/BMCodeAcc.h>
#include <fstream>

using namespace Rcpp;
using namespace std;

/******************************************************************************/

// [[Rcpp::export]]
bool readbina(const char * filename,
              Environment BM,
              const arma::Mat<unsigned char>& tab) {

  XPtr<FBM> xpBM = BM["address"];
  unsigned char* ptr = static_cast<unsigned char*>(xpBM->matrix());
  const unsigned char* code_ptr;
  int n = xpBM->nrow();
  int m = xpBM->ncol();

  int length = n / 4;
  int extra = n - 4 * length;
  int lengthExtra = length + (extra > 0);
  int j, k;

  unsigned char *buffer = new unsigned char[std::max(3, lengthExtra)];

  ifstream myFile(filename, ios::in | ios::binary);
  // magic number
  myFile.read((char*)buffer, 3);
  myassert((buffer[0] == 108) && (buffer[1] == 27) && (buffer[2] = 1),
           "Wrong magic number. Aborting..");

  for (j = 0; j < m; j++) {
    // read from bedfile
    myFile.read((char*)buffer, lengthExtra);

    for (k = 0; k < length; k++) {
      code_ptr = tab.colptr(buffer[k]);
      ptr = std::copy(code_ptr, code_ptr + 4, ptr);
    }
    if (extra) {
      code_ptr = tab.colptr(buffer[k]);
      ptr = std::copy(code_ptr, code_ptr + extra, ptr);
    }
  }

  char c;
  bool is_eof = !(myFile.get(c));
  myFile.close();
  delete[] buffer;

  return is_eof;
}

/******************************************************************************/

// [[Rcpp::export]]
void writebina(const char * filename,
               Environment BM,
               const RawVector& tab,
               const IntegerVector& rowInd,
               const IntegerVector& colInd) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd - 1, colInd - 1, BM["code256"]);
  int n = macc.nrow();
  int m = macc.ncol();
  int length = ceil((double)n / 4); // DO NOT USE INTEGERS WITH CEIL

  char *buffer = new char[std::max(3, length)];
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
