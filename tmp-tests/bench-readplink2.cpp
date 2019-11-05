#include <bigstatsr/BMAcc.h>
#include "../src/bed-acc.h"
#include <fstream>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
void fill_FBM(Environment BM,
              Environment obj_bed,
              const IntegerVector& ind_row,
              const IntegerVector& ind_col) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAcc macc_bed(xp_bed, ind_row, ind_col);

  XPtr<FBM> xpBM = BM["address"];
  BMAcc<unsigned char> macc_fbm(xpBM);

  size_t n = macc_bed.nrow();
  size_t m = macc_bed.ncol();

  for (size_t j = 0; j < m; j++)
    for (size_t i = 0; i < n; i++)
      macc_fbm(i, j) = macc_bed(i, j);
}

// [[Rcpp::export]]
bool fill_FBM2(const char * filename,
               Environment BM,
               const arma::Mat<unsigned char>& tab) {

  XPtr<FBM> xpBM = BM["address"];
  BMAcc<unsigned char> macc_fbm(xpBM);
  int n = xpBM->nrow();
  int m = xpBM->ncol();

  int length = n / 4;
  int extra = n - 4 * length;
  int lengthExtra = length + (extra > 0);
  int i, j, k;

  unsigned char *buffer = new unsigned char[std::max(3, lengthExtra) + 1];
  unsigned char *geno   = new unsigned char[4 * lengthExtra];
  const unsigned char *code_ptr;

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
      geno = std::copy(code_ptr, code_ptr + 4, geno);
    }
    if (extra) {
      code_ptr = tab.colptr(buffer[k]);
      geno = std::copy(code_ptr, code_ptr + extra, geno);
    }
    geno -= n;  // replace at start

    for (i = 0; i < n; i++) macc_fbm(i, j) = geno[i];
  }

  char c;
  bool is_eof = !(myFile.get(c));
  myFile.close();
  delete[] buffer;
  delete[] geno;

  return is_eof;
}

// [[Rcpp::export]]
void fill_FBM3(Environment BM,
               Environment obj_bed,
               const arma::Mat<unsigned char>& tab,
               ListOf<IntegerVector> ind) {

  XPtr<bed> xp_bed = obj_bed["address"];

  XPtr<FBM> xpBM = BM["address"];
  BMAcc<unsigned char> macc_fbm(xpBM);
  int n = xpBM->nrow();
  size_t m = xpBM->ncol();

  int length = n / 4;
  int extra = n - 4 * length;
  int k;
  size_t nbyte = xp_bed->nbyte();

  const unsigned char *code_ptr;
  const unsigned char *pMat0 = xp_bed->matrix();
  const unsigned char *pMat;

  // magic number
  myassert((pMat0[0] == 108) && (pMat0[1] == 27) && (pMat0[2] == 1),
           "Wrong magic number. Aborting..");
  pMat0 += 3;

  for (size_t j = 0; j < m; j++) {

    if (j % 1000 == 0) Rcout << j << " ";

    pMat = pMat0 + j * nbyte;

    for (k = 0; k < length; k++) {
      code_ptr = tab.colptr(pMat[k]);
      for (int l = 0; l < 3; l++) {
        IntegerVector ind_k = ind[k];
        if (ind_k.size()) {
          const unsigned char geno = code_ptr[l];
          for (auto& i : ind_k) macc_fbm(i, j) = geno;
        }
      }
    }

    if (extra) {
      code_ptr = tab.colptr(pMat[k]);
      for (int l = 0; l < extra; l++) {
        IntegerVector ind_k = ind[k];
        if (ind_k.size()) {
          const unsigned char geno = code_ptr[l];
          for (auto& i : ind_k) macc_fbm(i, j) = geno;
        }
      }
    }
  }
}

// [[Rcpp::export]]
void fill_FBM4(Environment BM,
               Environment obj_bed,
               const RawMatrix& tab,
               const IntegerVector& size,
               const IntegerVector& ind) {

  XPtr<bed> xp_bed = obj_bed["address"];

  XPtr<FBM> xpBM = BM["address"];
  BMAcc<unsigned char> macc_fbm(xpBM);
  int n = xpBM->nrow();
  size_t m = xpBM->ncol();

  int length = n / 4;
  int extra = n - 4 * length;
  int k;
  size_t nbyte = xp_bed->nbyte();

  const unsigned char *pMat0 = xp_bed->matrix();
  const unsigned char *pMat;

  // magic number
  myassert((pMat0[0] == 108) && (pMat0[1] == 27) && (pMat0[2] == 1),
           "Wrong magic number. Aborting..");
  pMat0 += 3;

  for (size_t j = 0; j < m; j++) {

    // if (j % 10000 == 0) Rcout << j << " ";

    pMat = pMat0 + j * nbyte;
    size_t i_size = 0, i_ind = 0;

    for (k = 0; k < length; k++) {
      int k2 = pMat[k];
      for (int l = 0; l < 3; l++) {
        int size_i = size[i_size++];
        if (size_i) {
          const unsigned char geno = tab(l, k2);
          for (int i = 0; i < size_i; i++) macc_fbm(ind[i_ind++], j) = geno;
        }
      }
    }

    if (extra) {
      int k2 = pMat[k];
      for (int l = 0; l < extra; l++) {
        int size_i = size[i_size++];
        if (size_i) {
          const unsigned char geno = tab(l, k2);
          for (int i = 0; i < size_i; i++) macc_fbm(ind[i_ind++], j) = geno;
        }
      }
    }
  }
}

/*** R
library(bigsnpr)

dim(obj.bed <- bed("tmp-data/1000G_phase3_common_hapmap_norel.bed"))
system.time(snp_readBed(obj.bed$bedfile, "tmp-data/test1"))
# 10.980   3.701  34.891
# 11.296   3.885  49.410
# 8.510   3.043  37.778

system.time(snp_readBed2(obj.bed$bedfile, "tmp-data/test2"))

G2 <- FBM(nrow(obj.bed), ncol(obj.bed), type = "raw",
          backingfile = "tmp-data/test6")
system.time({
  fill_FBM(G2, obj.bed, rows_along(obj.bed), cols_along(obj.bed))
  G2[, 1]
})
# 11.226   3.833  31.821
# 9.158   3.157  31.468
# 8.206   3.285  34.881

G <- FBM(nrow(obj.bed), ncol(obj.bed), type = "raw",
         backingfile = "tmp-data/test4")
system.time(fill_FBM2(obj.bed$bedfile, G, bigsnpr:::getCode()))
# 8.135   3.674  35.118

G3 <- FBM(nrow(obj.bed), ncol(obj.bed), type = "raw",
          backingfile = "tmp-data/test10")
# system.time(fill_FBM3(G3, obj.bed, bigsnpr:::getCode(),
#                       as.list(rows_along(obj.bed) - 1)))
system.time(fill_FBM4(G3, obj.bed, bigsnpr:::getCode(),
                      rep(1, nrow(G3)), rows_along(G3) - 1))
# 7.021   3.008  53.181
*/
