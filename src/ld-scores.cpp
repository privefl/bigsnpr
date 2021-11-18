/******************************************************************************/

// #include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/BMCodeAcc.h>
#include "bed-acc.h"

using namespace Rcpp;

/******************************************************************************/

template <class C>
NumericVector ld_scores0(C macc,
                         double size,
                         const NumericVector& pos,
                         int ncores) {

  int n = macc.nrow();
  int m = macc.ncol();

  std::vector<double> res(m, 1);

  int chunk_size = ceil(m / (10.0 * ncores));

  #pragma omp parallel for schedule(dynamic, chunk_size) num_threads(ncores)
  for (int j0 = 0; j0 < m; j0++) {

    // pre-computation
    double xSum0 = 0, xxSum0 = 0;
    for (int i = 0; i < n; i++) {
      double x = macc(i, j0);
      if (x != 3) {
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
        if (x == 3) continue;

        double y = macc(i, j);
        if (y == 3) {
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
      double r2 = num * num / (deno_x * deno_y);

      if (!ISNAN(r2)) {
        #pragma omp atomic
        res[j0] += r2;
        #pragma omp atomic
        res[j] += r2;
      }
    }
  }

  return Rcpp::wrap(res);
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector ld_scores(Environment obj,
                        const IntegerVector& rowInd,
                        const IntegerVector& colInd,
                        double size,
                        const NumericVector& pos,
                        int ncores) {

  myassert_size(colInd.size(), pos.size());

  if (obj.exists("code256")) {
    XPtr<FBM> xpBM = obj["address"];
    NumericVector code = clone(as<NumericVector>(obj["code256"]));
    code[is_na(code)] = 3;
    SubBMCode256Acc macc(xpBM, rowInd, colInd, code, 1);
    return ld_scores0(macc, size, pos, ncores);
  } else if (obj.exists("bedfile")) {
    XPtr<bed> xp_bed = obj["address"];
    bedAcc macc(xp_bed, rowInd, colInd);
    return ld_scores0(macc, size, pos, ncores);
  } else {
    throw Rcpp::exception("Unknown object type.");
  }
}

/******************************************************************************/
