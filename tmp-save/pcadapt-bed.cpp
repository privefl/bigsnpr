/******************************************************************************/

#include "bed-acc.hpp"

/******************************************************************************/

// [[Rcpp::export]]
NumericMatrix multLinReg(Environment obj_bed,
                         const IntegerVector& ind_row,
                         const IntegerVector& ind_col,
                         const NumericVector& center,
                         const NumericVector& scale,
                         const NumericMatrix& u) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAccScaled macc(xp_bed, ind_row, ind_col, center, scale);
  size_t n = macc.nrow();
  size_t p = macc.ncol();
  int K = u.ncol();

  NumericMatrix Z(p, K);
  NumericMatrix x(n);
  double eps;
  int k, nona;

  for (size_t j = 0; j < p; j++) {

    LogicalVector not_missing(n);

    // counting NAs and Z = U^T G (x)
    nona = n;
    NumericVector z(K);              // all 0s
    for (size_t i = 0; i < n; i++) {
      x[i] = macc(i, j);
      not_missing[i] = (x[i] != 3);
      if (not_missing[i]) {
        for (k = 0; k < K; k++) {
          z[k] += u(i, k) * x[i];
        }
      } else {
        nona--;
      }
    }

    // G* (y) = U Z
    NumericVector y(n);              // all 0s
    NumericVector sum_scores_sq(K);  // all 0s
    double sum_resid_sq = 0;
    for (size_t i = 0; i < n; i++) {
      if (not_missing[i]) {
        for (k = 0; k < K; k++) {
          y[i] += u(i, k) * z[k];
          sum_scores_sq[k] += u(i, k) * u(i, k);  // can't precompute
        }
        eps = x[i] - y[i];
        sum_resid_sq += eps * eps;
      }
    }

    for (k = 0; k < K; k++) {
      Z(j, k) = z[k] / sqrt(sum_scores_sq[k] * sum_resid_sq / (nona - K));
    }
  }

  return Z;
}

/******************************************************************************/
