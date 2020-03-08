/******************************************************************************/

#include "bed-acc.h"


#if defined(_OPENMP)
#include <omp.h>
#else
inline int omp_get_thread_num() { return 0; }
#endif

/******************************************************************************/

// [[Rcpp::export]]
NumericVector bed_pMatVec4(Environment obj_bed,
                           const IntegerVector& ind_row,
                           const IntegerVector& ind_col,
                           const NumericVector& center,
                           const NumericVector& scale,
                           const NumericVector& x,
                           int ncores) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAccScaled macc(xp_bed, ind_row, ind_col, center, scale);

  int n = macc.nrow();
  NumericMatrix res(n, ncores);

  #pragma omp parallel num_threads(ncores)
  {
    int id = omp_get_thread_num();
    int n2 = n;
    int m = macc.ncol();
    int m2 = m - 3;  // WARNING: do not use std::size_t because of `m - 3`
    int i, j;

    // unrolling optimization
    #pragma omp for nowait
    for (j = 0; j < m2; j += 4) {
      for (i = 0; i < n2; i++) {
        res(i, id) += (x[j] * macc(i, j) + x[j+1] * macc(i, j+1)) +
          (x[j+2] * macc(i, j+2) + x[j+3] * macc(i, j+3));
      } // The parentheses are somehow important
    }
    #pragma omp for
    for (j = m - m % 4; j < m; j++) {
      for (i = 0; i < n2; i++) {
        res(i, id) += x[j] * macc(i, j);
      }
    }
  }

  return rowSums(res);
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector bed_cpMatVec4(Environment obj_bed,
                            const IntegerVector& ind_row,
                            const IntegerVector& ind_col,
                            const NumericVector& center,
                            const NumericVector& scale,
                            const NumericVector& x,
                            int ncores) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAccScaled macc(xp_bed, ind_row, ind_col, center, scale);

  int m = macc.ncol();
  NumericVector res(m);

  #pragma omp parallel num_threads(ncores)
  {
    int m2 = m;
    int n = macc.nrow();
    int n2 = n - 3;  // WARNING: do not use std::size_t because of `n - 3`

    #pragma omp for
    for (int j = 0; j < m2; j++) {

      double tmp = 0;

      // unrolling optimization
      int i = 0;
      for (; i < n2; i += 4) {
        tmp += (macc(i, j) * x[i] + macc(i+1, j) * x[i+1]) +
          (macc(i+2, j) * x[i+2] + macc(i+3, j) * x[i+3]);
      }
      for (; i < n; i++) tmp += macc(i, j) * x[i];

      res[j] = tmp;
    }
  }

  return res;
}

/******************************************************************************/
