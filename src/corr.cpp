/******************************************************************************/

// #include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/BMCodeAcc.h>
#include "bed-acc.h"

using namespace Rcpp;

/******************************************************************************/

template <class C>
List corMat0(C macc,
             double size,
             const NumericVector& thr,
             const NumericVector& pos,
             bool fill_diag,
             int ncores) {

  int n = macc.nrow();
  int m = macc.ncol();

  List res(m);

  int chunk_size = ceil(m / (10.0 * ncores));

  #pragma omp parallel num_threads(ncores)
  {
    std::vector<int>    ind; ind.reserve(m);
    std::vector<double> val; val.reserve(m);

    #pragma omp for schedule(dynamic, chunk_size)
    for (int j0 = 0; j0 < m; j0++) {

      ind.clear();
      val.clear();
      if (fill_diag) {
        ind.push_back(j0);
        val.push_back(1.0);
      }

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
        double r = num / ::sqrt(deno_x * deno_y);

        if (ISNAN(r) || std::abs(r) > thr[nona - 1]) {
          ind.push_back(j);
          val.push_back(r);
        }
      }

      #pragma omp critical
      res[j0] = List::create(
        _["i"] = rev(as<IntegerVector>(wrap(ind))),
        _["x"] = rev(as<NumericVector>(wrap(val))));
    }
  }

  return res;
}

/******************************************************************************/

// [[Rcpp::export]]
List corMat(Environment obj,
            const IntegerVector& rowInd,
            const IntegerVector& colInd,
            double size,
            const NumericVector& thr,
            const NumericVector& pos,
            bool fill_diag,
            int ncores) {

  myassert_size(colInd.size(), pos.size());

  if (obj.exists("code256")) {
    XPtr<FBM> xpBM = obj["address"];
    NumericVector code = clone(as<NumericVector>(obj["code256"]));
    code[is_na(code)] = 3;
    SubBMCode256Acc macc(xpBM, rowInd, colInd, code, 1);
    return corMat0(macc, size, thr, pos, fill_diag, ncores);
  } else if (obj.exists("bedfile")) {
    XPtr<bed> xp_bed = obj["address"];
    bedAcc macc(xp_bed, rowInd, colInd);
    return corMat0(macc, size, thr, pos, fill_diag, ncores);
  } else {
    throw Rcpp::exception("Unknown object type.");
  }
}

/******************************************************************************/
