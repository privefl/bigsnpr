/******************************************************************************/

#include "bigsnpr.h"
#include <bigmemory/isna.hpp>

/******************************************************************************/

// [[Rcpp::export]]
SEXP corMat(const S4& BM,
            const IntegerVector& rowInd,
            const IntegerVector& colInd,
            int size,
            const NumericVector& thr) {

  XPtr<BigMatrix> xpMat = BM.slot("address");
  RawSubMatAcc macc(*xpMat, rowInd-1, colInd-1, BM.slot("code"));

  int n = macc.nrow();
  int m = macc.ncol();

  arma::sp_mat corr(m, m);

  int i, j, j0, N;
  double x, y;
  double xSum, xxSum, deno_x;
  double ySum, yySum, deno_y;
  double xySum, num, r, t;

  // pre-computation
  NumericVector sumX(m), sumXX(m);

  for (j = 0; j < m; j++) {
    xSum = xxSum = 0;
    for (i = 0; i < n; i++) {
      x = macc(i, j);

      if (!isna(x)) {
        xSum += x;
        xxSum += x*x;
      }
    }
    sumX[j] = xSum;
    sumXX[j] = xxSum;
  }

  // main computation
  for (j0 = 0; j0 < m; j0++) {
    for (j = max(0, j0 - size); j < j0; j++) {
      N = n;
      xySum = 0;
      xSum = sumX[j];
      ySum = sumX[j0];
      xxSum = sumXX[j];
      yySum = sumXX[j0];
      for (i = 0; i < n; i++) {
        x = macc(i, j);
        y = macc(i, j0);

        if (isna(y)) {
          // printf("Missing value at pos (%d, %d)\n", i+1, j0+1); // DEBUG
          N--;
          if (isna(x)) { // both missing
            // nothing to do
          } else { // y is missing but not x
            xSum -= x;
            xxSum -= x*x;
          }
        } else {
          if (isna(x)) { // x is missing but not y
            ySum -= y;
            yySum -= y*y;
            N--;
          } else { // both not missing
            xySum += x * y;
          }
        }
      }
      // printf("N = %d\n", N); // DEBUG
      num = xySum - xSum * ySum / N;
      deno_x = xxSum - xSum * xSum / N;
      deno_y = yySum - ySum * ySum / N;
      r = num / sqrt(deno_x * deno_y);
      t = r * sqrt((N - 2) / (1 - r*r));
      // printf("%f -- %d -- %f -- %f\n", r, N, t, thr[N-1]); // DEBUG
      if (abs(t) > thr[N-1]) corr(j, j0) = r;
    }
  }

  return wrap(corr);
}

/******************************************************************************/
