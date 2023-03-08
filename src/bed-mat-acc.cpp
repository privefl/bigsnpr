/******************************************************************************/

#include "bed-acc.h"

/******************************************************************************/

// [[Rcpp::export]]
IntegerMatrix read_bed(Environment obj_bed,
                       const IntegerVector& ind_row,
                       const IntegerVector& ind_col) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAcc macc_bed(xp_bed, ind_row, ind_col, NA_INTEGER);

  size_t n = macc_bed.nrow();
  size_t m = macc_bed.ncol();

  IntegerMatrix res(n, m);

  for (size_t j = 0; j < m; j++)
    for (size_t i = 0; i < n; i++)
      res(i, j) = macc_bed(i, j);

  return res;
}

/******************************************************************************/

// [[Rcpp::export]]
NumericMatrix read_bed_scaled(Environment obj_bed,
                              const IntegerVector& ind_row,
                              const IntegerVector& ind_col,
                              const NumericVector& center,
                              const NumericVector& scale) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAccScaled macc_bed(xp_bed, ind_row, ind_col, center, scale);

  size_t n = macc_bed.nrow();
  size_t m = macc_bed.ncol();

  NumericMatrix res(n, m);

  for (size_t j = 0; j < m; j++)
    for (size_t i = 0; i < n; i++)
      res(i, j) = macc_bed(i, j);

  return res;
}

/******************************************************************************/
