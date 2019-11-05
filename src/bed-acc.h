#ifndef BED_ACC_H
#define BED_ACC_H

/******************************************************************************/

#include <bigstatsr/utils.h>

using namespace Rcpp;
using std::size_t;

/******************************************************************************/

class bed {
public:
  bed(const std::string path, int n, int m);

  IntegerMatrix get_code(int NA_VAL = 3) const {

    IntegerVector num = IntegerVector::create(2, NA_VAL, 1, 0);
    IntegerMatrix code(4, 256);

    int i, k, k2;
    int coeff = 1;
    for (i = 0; i < 4; i++) {
      for (k = 0; k < 256; k++) {
        k2 = k / coeff;
        code(i, k) = num[k2 % 4];
      }
      coeff *= 4;
    }

    return code;
  }

  const unsigned char* matrix() const { return ro_ummap.data(); }

  size_t ntot() const { return n; }
  size_t mtot() const { return m; }
  size_t nbyte() const { return n_byte; }

private:
  mio::ummap_source ro_ummap;
  size_t n, m, n_byte;
};

/******************************************************************************/

class bedAcc {
public:
  bedAcc(bed * bedPtr,
         const IntegerVector& ind_row,
         const IntegerVector& ind_col) {

    size_t n = bedPtr->ntot();
    size_t m = bedPtr->mtot();
    n_byte = bedPtr->nbyte();
    _pMat = 3 + bedPtr->matrix();
    _lookup_byte = bedPtr->get_code();

    // Indices of sub-view of bed file
    size_t n_sub = ind_row.size();
    std::vector<size_t> ind_row2(n_sub);
    for (size_t i = 0; i < n_sub; i++) {
      // 'ind_row' indices comes from R, so begins at 1
      size_t ind = static_cast<size_t>(ind_row[i] - 1);
      myassert(ind < n, ERROR_BOUNDS);
      ind_row2[i] = ind;
    }
    _ind_row = ind_row2;

    size_t m_sub = ind_col.size();
    std::vector<size_t> ind_col2(m_sub);
    for (size_t j = 0; j < m_sub; j++) {
      // 'ind_col' indices comes from R, so begins at 1
      size_t ind = static_cast<size_t>(ind_col[j] - 1);
      myassert(ind < m, ERROR_BOUNDS);
      ind_col2[j] = ind;
    }
    _ind_col = ind_col2;
  }

  size_t nrow() const { return _ind_row.size(); }
  size_t ncol() const { return _ind_col.size(); }

  inline int operator() (size_t i, size_t j) const {
    size_t i2 = _ind_row[i];
    unsigned char byte = _pMat[i2 / 4 + _ind_col[j] * n_byte];
    return _lookup_byte(i2 % 4, byte);
  }

protected:
  size_t n_byte;
  const unsigned char* _pMat;
  IntegerMatrix _lookup_byte;
  std::vector<size_t> _ind_row, _ind_col;
};

/******************************************************************************/

class bedAccScaled : public bedAcc {
public:
  bedAccScaled(bed * bedPtr,
               const IntegerVector& ind_row,
               const IntegerVector& ind_col,
               const NumericVector& center,
               const NumericVector& scale,
               double NA_VAL = 0) : bedAcc(bedPtr, ind_row, ind_col) {

    myassert_size(center.size(), ind_col.size());
    myassert_size(scale.size(),  ind_col.size());

    size_t p = ind_col.size();
    _lookup_scale = NumericMatrix(4, p);
    for (size_t j = 0; j < p; j++) {
      for (size_t i = 0; i < 3; i++) {
        _lookup_scale(i, j) = (i - center[j]) / scale[j];
      }
      _lookup_scale(3, j) = NA_VAL;
    }
  };

  inline double operator() (size_t i, size_t j) const {
    int geno = bedAcc::operator()(i, j);
    return _lookup_scale(geno, j);
  }

protected:
  NumericMatrix _lookup_scale;
};

/******************************************************************************/

#endif // BED_ACC_H
