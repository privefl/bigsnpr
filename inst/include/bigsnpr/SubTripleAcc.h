#ifndef SUBTRIPLEACC_H
#define SUBTRIPLEACC_H

/******************************************************************************/

#include <bigstatsr/SubMatAcc.h>

/******************************************************************************/

class RawSubTripleAcc : public SubMatAcc<unsigned char> {
public:
  RawSubTripleAcc(BigMatrix& bm,
                  const IntegerVector& row_ind,
                  const IntegerVector& col_ind,
                  const NumericVector& lookup)
    : SubMatAcc<unsigned char>(bm, row_ind, col_ind) {
      NumericMatrix tmp(lookup.size(), 3);
      tmp(_, 0) = lookup;
      tmp(_, 1) = lookup >= 0.5;
      tmp(_, 2) = lookup > 1.5;
      _lookup = tmp;
    }

  inline double operator() (int i, int j) {
    int j2 = j / 3;
    return _lookup(*(_pMat + _totalRows * _col_ind[j2] + _row_ind[i]), j % 3);
  }

  int ncol() const {
    return _ncol * 3;
  }

protected:
  NumericMatrix _lookup;
};

/******************************************************************************/

#endif // SUBTRIPLEACC_H
