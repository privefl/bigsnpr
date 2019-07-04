#ifndef BED_ACC_H
#define BED_ACC_H

/******************************************************************************/

#include <bigstatsr/utils.h>
#include <system_error> // for std::error_code

using namespace Rcpp;
using std::size_t;

/******************************************************************************/

class bedAcc {
public:
  bedAcc(const std::string path, size_t n, size_t m,
         const IntegerVector& row_ind,
         const IntegerVector& col_ind,
         const RawMatrix& lookup_byte) :
  n(n), m(m), n_byte((n + 3) / 4), _lookup_byte(lookup_byte) {

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
    if ((this->n_byte * this->m) != this->ro_ummap.size())
      Rcpp::stop("n or m does not match the dimensions of the file.");

    _pMat = ro_ummap.data();


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
      size_t ind = static_cast<size_t>(col_ind[j] - 1);
      myassert(ind < m, ERROR_BOUNDS);
      col_ind2[j] = ind;
    }
    _col_ind = col_ind2;
  };

  size_t nrow() const { return _row_ind.size(); }
  size_t ncol() const { return _col_ind.size(); }

  inline unsigned char operator() (size_t i, size_t j) {
    size_t i2 = _row_ind[i];
    const unsigned char byte = _pMat[i2 / 4 + _col_ind[j] * n_byte];
    return _lookup_byte(i2 % 4, byte);
  }

protected:
  size_t n, m, n_byte;
  mio::ummap_source ro_ummap;
  const unsigned char* _pMat;
  RawMatrix _lookup_byte;
  std::vector<size_t> _row_ind, _col_ind;
};

/******************************************************************************/

class bedAccScaled : public bedAcc {
public:
  bedAccScaled(const std::string path, size_t n, size_t m,
               const IntegerVector& row_ind,
               const IntegerVector& col_ind,
               const RawMatrix& lookup_byte,
               const NumericMatrix& lookup_scale) :
  bedAcc(path, n, m, row_ind, col_ind, lookup_byte) {

    myassert_size(lookup_scale.nrow(), 4);
    myassert_size(lookup_scale.ncol(), col_ind.size());

    _lookup_scale = lookup_scale;
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
