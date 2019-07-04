#ifndef BED_ACC_H
#define BED_ACC_H

/******************************************************************************/

#include <mio/mmap.hpp>
#include <system_error> // for std::error_code
#include <Rcpp.h>

using namespace Rcpp;
using std::size_t;

/******************************************************************************/

class bedAcc {
public:
  bedAcc(const std::string path, int n, int p,
         const IntegerVector& row_ind,
         const IntegerVector& col_ind,
         const RawMatrix& lookup_byte) {

    size_t n = n, p = p, n_byte = (n + 3) / 4;
    _lookup_byte = lookup_byte;


    // Memory-map the bed file
    std::error_code error;
    mio::ummap_source magic_number;
    magic_number.map(path, 0, 3, error);
    if (error) Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

    if (!(magic_number[0] == '\x6C' && magic_number[1] == '\x1B'))
      Rcpp::stop("File is not a bed file.");

    /* Check mode: 00000001 indicates the default variant-major mode (i.e.
    list all samples for first variant, all samples for second variant,
    etc), 00000000 indicates the unsupported sample-major mode (i.e. list
    all variants for the first sample, list all variants for the second
    sample, etc */
    if (magic_number[2] != '\x01') Rcpp::stop("Sample-major mode is not supported.");

    // Map after this magic number
    this->ro_ummap.map(path, 3, mio::map_entire_file, error);
    if (error) Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

    // Check if given dimensions match the file
    if ((this->n_byte * this->p) != this->ro_ummap.size())
      Rcpp::stop("n or p does not match the dimensions of the file.");


    // Indices of sub-view of bed file
    size_t n_sub = row_ind.size();
    std::vector<size_t> row_ind2(n_sub);
    for (size_t i = 0; i < n_sub; i++) {
      // 'row_ind' indices comes from R, so begins at 1
      size_t ind = static_cast<size_t>(row_ind[i] - 1);
      myassert(ind < n, ERROR_BOUNDS);
      row_ind2[i] = ind;
    }
    _row_ind = row_ind2;

    size_t p_sub = col_ind.size();
    std::vector<size_t> col_ind2(p_sub);
    for (size_t j = 0; j < p_sub; j++) {
      // 'col_ind' indices comes from R, so begins at 1
      ind = static_cast<size_t>(col_ind[j] - 1);
      myassert(ind < p, ERROR_BOUNDS);
      col_ind2[j] = ind;
    }
    _col_ind = col_ind2;
  };

  size_t nrow() const { return _row_ind.size(); }
  size_t ncol() const { return _col_ind.size(); }

  inline unsigned char operator() (size_t i, size_t j) {
    int i2 = _row_ind[i];
    const unsigned char byte = _pMat[i / 4 + _col_ind[j] * n_byte];
    return _lookup_byte(i % 4, byte);
  }

protected:
  size_t n, p, n_byte;
  mio::ummap_source ro_ummap;
  const unsigned char* _pMat;
  RawMatrix _lookup_byte;
  std::vector<size_t> _row_ind, _col_ind;
};

/******************************************************************************/

class bedAccScaled : public bedAcc {
public:
  bedAccScaled(const std::string path, int n, int p,
               const IntegerVector& row_ind,
               const IntegerVector& col_ind,
               const IntegerVector& lookup_byte,
               // af should be ALL allele frequencies
               const NumericVector& af,
               double NA_VAL,
               double ploidy = 2) :
  bedAcc(path, n, p, row_ind, col_ind, lookup_byte) {

    _lookup_scale = NumericMatrix(4, p);
    for (size_t j = 0; j < p; j++) {
      double af_j = af[_col_ind[j]];
      for (size_t i = 0; i < 3; i++) {
        _lookup_scale(i, j) =
          (i - ploidy * af_j) / sqrt(ploidy * af_j * (1 - af_j));
      }
      _lookup_scale(3, j) = NA_VAL;
    }
  };

  inline double operator() (size_t i, size_t j) {
    unsigned char geno = bedAcc::operator()(i, j);
    return _lookup_scale(geno, j);
  }

protected:
  NumericMatrix _lookup_scale;
};

/******************************************************************************/

#endif // BED_ACC_H
