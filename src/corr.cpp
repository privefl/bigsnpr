/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
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
            const NumericVector& pos,
            int ncores = 1) {

  myassert_size(colInd.size(), pos.size());

  XPtr<FBM> xpBM = BM["address"];
  NumericVector code = clone(as<NumericVector>(BM["code256"]));
  code[is_na(code)] = 3;
  SubBMCode256Acc macc(xpBM, rowInd, colInd, code, 1);

  int n = macc.nrow();
  int m = macc.ncol();

  arma::sp_mat corr(m, m);

  #pragma omp parallel for num_threads(ncores)
  for (int j0 = 0; j0 < m; j0++) {

    // pre-computation
    double xSum0 = 0, xxSum0 = 0;
    for (int i = 0; i < n; i++) {
      double x = macc(i, j0);
      if (!isna(x)) {
        xSum0  += x;
        xxSum0 += x * x;
      }
    }

    // main computation
    double pos_min = pos[j0] - size;
    for (int j = j0 - 1; (j >= 0) && (pos[j] >= pos_min); j--) {

      int nona = 0;
      double xSum = xSum0, xxSum = xxSum0;
      double ySum = 0, yySum = 0, xySum = 0;
      for (int i = 0; i < n; i++) {

        double x = macc(i, j0);
        if (isna(x)) continue;

        double y = macc(i, j);
        if (isna(y)) {
          // y is missing, but not x
          xSum  -= x;
          xxSum -= x * x;
        } else {
          // none missing
          nona++;
          ySum  += y;
          yySum += y * y;
          xySum += x * y;
        }
      }

      double num = xySum - xSum * ySum / nona;
      double deno_x = xxSum - xSum * xSum / nona;
      double deno_y = yySum - ySum * ySum / nona;
      double r = num / ::sqrt(deno_x * deno_y);

      if (std::abs(r) > thr[nona - 1]) {
        #pragma omp critical
        corr(j, j0) = r;
      }
    }
  }

  return wrap(corr);
}

/******************************************************************************/
