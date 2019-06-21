/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

using namespace Rcpp;

/******************************************************************************/

// R_IsNA:  https://stackoverflow.com/a/26262984/6103040
// Using 3: https://stackoverflow.com/q/46892399/6103040
inline bool isna(double x) {
  return(x == 3);
}

// [[Rcpp::export]]
SEXP corMat(Environment BM,
            const IntegerVector& rowInd,
            const IntegerVector& colInd,
            double size,
            const NumericVector& thr,
            const NumericVector& pos) {

  XPtr<FBM> xpBM = BM["address"];
  NumericVector code = clone(as<NumericVector>(BM["code256"]));
  code[is_na(code)] = 3;
  SubBMCode256Acc macc(xpBM, rowInd - 1, colInd - 1, code);

  int n = macc.nrow();
  int m = macc.ncol();

  arma::sp_mat corr(m, m);

  int i, j, j0, N;
  double x, y;
  double xSum, xxSum, deno_x;
  double ySum, yySum, deno_y;
  double pos_min, xySum, num, r;

  NumericVector sumX(m), sumXX(m);

  for (j0 = 0; j0 < m; j0++) {

    // pre-computation
    xSum = xxSum = 0;
    for (i = 0; i < n; i++) {
      x = macc(i, j0);
      if (!isna(x)) {
        xSum += x;
        xxSum += x*x;
      }
    }
    sumX[j0] = xSum;
    sumXX[j0] = xxSum;

    // main computation
    pos_min = pos[j0] - size;
    for (j = j0 - 1; (j >= 0) && (pos[j] >= pos_min); j--) {

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
          N--;
          if (isna(x)) {       // both missing
            // nothing to do
          } else {             // y is missing, but not x
            xSum -= x;
            xxSum -= x*x;
          }
        } else {
          if (isna(x)) {       // x is missing, but not y
            ySum -= y;
            yySum -= y*y;
            N--;
          } else {             // none missing
            xySum += x * y;
          }
        }
      }

      num = xySum - xSum * ySum / N;
      deno_x = xxSum - xSum * xSum / N;
      deno_y = yySum - ySum * ySum / N;
      r = num / sqrt(deno_x * deno_y);

      if (std::abs(r) > thr[N-1]) corr(j, j0) = r;
    }
  }

  return wrap(corr);
}

/******************************************************************************/
