/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

using namespace Rcpp;

/******************************************************************************/

bool is_na0(double x) { return NumericVector::is_na(x); }

bool is_na3(double x) { return(x == 3); }

/******************************************************************************/

template <class C>
List corMat0(C macc,
             bool(*isna)(double),
             const IntegerVector& rowInd,
             const IntegerVector& colInd,
             double size,
             const NumericVector& thr,
             const NumericVector& pos,
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
          ind.push_back(j + 1);
          val.push_back(r);
        }
      }

      #pragma omp critical
      res[j0] = List::create(_["i"] = wrap(ind), _["x"] = wrap(val));
    }
  }


  return res;
}

/******************************************************************************/

// [[Rcpp::export]]
List corMat(Environment BM,
            const IntegerVector& rowInd,
            const IntegerVector& colInd,
            double size,
            const NumericVector& thr,
            const NumericVector& pos,
            int ncores) {

  myassert_size(colInd.size(), pos.size());

  XPtr<FBM> xpBM = BM["address"];

  if (BM.exists("code256")) {
    NumericVector code = clone(as<NumericVector>(BM["code256"]));
    code[is_na(code)] = 3;
    SubBMCode256Acc macc(xpBM, rowInd, colInd, code, 1);
    return corMat0(macc, &is_na3, rowInd, colInd, size, thr, pos, ncores);
  } else {
    switch(xpBM->matrix_type()) {
    case 6:
    {
      SubBMAcc<float> macc(xpBM, rowInd, colInd, 1);
      return corMat0(macc, &is_na0, rowInd, colInd, size, thr, pos, ncores);
    }
    default:
      throw Rcpp::exception(ERROR_TYPE);
    }
  }
}

/******************************************************************************/
