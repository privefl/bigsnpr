/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <fstream>
#include "bed-acc.h"

using namespace Rcpp;
using namespace std;

/******************************************************************************/

// [[Rcpp::export]]
bool readbina(const char * filename,
              Environment BM,
              const RawMatrix& tab) {

  XPtr<FBM_RW> xpBM = BM["address_rw"];
  unsigned char* ptr = static_cast<unsigned char*>(xpBM->matrix());
  const unsigned char* code_ptr;
  int n = xpBM->nrow();
  int m = xpBM->ncol();

  int length = n / 4;
  int extra = n - 4 * length;
  int lengthExtra = length + (extra > 0);
  int j, k;

  unsigned char *buffer = new unsigned char[std::max(3, lengthExtra) + 1];

  ifstream myFile(filename, ios::in | ios::binary);
  // magic number
  myFile.read((char*)buffer, 3);
  myassert((buffer[0] == 108) && (buffer[1] == 27) && (buffer[2] = 1),
           "Wrong magic number. Aborting..");

  for (j = 0; j < m; j++) {
    // read from bedfile
    myFile.read((char*)buffer, lengthExtra);

    for (k = 0; k < length; k++) {
      code_ptr = &tab(0, buffer[k]);
      ptr = std::copy(code_ptr, code_ptr + 4, ptr);
    }
    if (extra) {
      code_ptr = &tab(0, buffer[k]);
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
void readbina2(Environment BM,
               Environment obj_bed,
               const IntegerVector& ind_row,
               const IntegerVector& ind_col,
               int ncores) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAcc macc_bed(xp_bed, ind_row, ind_col);

  XPtr<FBM_RW> xpBM = BM["address_rw"];
  BMAcc_RW<unsigned char> macc_fbm(xpBM);

  size_t n = macc_bed.nrow();
  size_t m = macc_bed.ncol();

  #pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++)
    for (size_t i = 0; i < n; i++)
      macc_fbm(i, j) = macc_bed(i, j);
}

/******************************************************************************/
