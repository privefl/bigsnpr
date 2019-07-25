/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include "bed-acc.h"

using namespace Rcpp;

/******************************************************************************/

template <class C>
NumericMatrix multLinReg(C macc, const NumericMatrix& U) {

  int n = macc.nrow();
  int m = macc.ncol();
  int K = U.ncol();
  myassert_size(U.nrow(), n);

  NumericMatrix res(m, K);

  for (int j = 0; j < m; j++) {
    NumericVector xySum(K), ySum(K), yySum(K);
    double xSum = 0, xxSum = 0;
    int nona = n;
    for (int i = 0; i < n; i++) {
      double x = macc(i, j);
      if (x != 3) { // not missing
        xSum  += x;
        xxSum += x*x;
        for (int k = 0; k < K; k++) {
          double y = U(i, k);
          xySum[k] += x * y;
          ySum[k]  += y;
          yySum[k] += y*y;
        }
      } else {
        nona--;
      }
    }

    double deno_x = xxSum - xSum * xSum / nona;
    for (int k = 0; k < K; k++) {
      double num = xySum[k] - xSum * ySum[k] / nona;
      double deno_y = yySum[k] - ySum[k] * ySum[k] / nona;
      double deno = deno_x * deno_y - num * num;
      res(j, k) = (deno == 0 || nona < 2) ? NA_REAL : num * sqrt((nona - 2) / deno);
    }
  }

  return res;  // t-scores
}

/******************************************************************************/

// Dispatch function for multLinReg
// [[Rcpp::export]]
NumericMatrix multLinReg(SEXP obj,
                         const IntegerVector& ind_row,
                         const IntegerVector& ind_col,
                         const NumericMatrix& U) {

  if (Rf_inherits(obj, "FBM.code256")) {
    Environment obj2 = obj;
    XPtr<FBM> xpBM = obj2["address"];
    NumericVector code = clone(as<NumericVector>(obj2["code256"]));
    code[is_na(code)] = 3;
    SubBMCode256Acc macc(xpBM, ind_row - 1, ind_col - 1, code);
    return multLinReg(macc, U);
  } else if (Rf_inherits(obj, "bed")) {
    XPtr<bed> xp_bed = as<Environment>(obj)["address"];
    bedAcc macc(xp_bed, ind_row, ind_col);
    return multLinReg(macc, U);
  } else {
    Rcpp::stop("multLinReg() is not implemented for this class.");
  }
}

/******************************************************************************/
