/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

std::vector<size_t>& which_to_check(int j0,
                                    const LogicalVector& keep,
                                    const NumericVector& pos,
                                    double size,
                                    std::vector<size_t>& ind_to_check) {

  ind_to_check.clear();
  int m = pos.size();

  int pos_min = pos[j0] - size;
  int pos_max = pos[j0] + size;
  bool not_min = true;
  bool not_max = true;

  for (int l = 1; not_max || not_min; l++) {
    if (not_max) {
      int j = j0 + l;
      not_max = (j < m) && (pos[j] <= pos_max);  // within a window..
      if (not_max && keep[j]) ind_to_check.push_back(j);
    }
    if (not_min) {
      int j = j0 - l;
      not_min = (j >= 0) && (pos[j] >= pos_min);  // within a window..
      if (not_min && keep[j]) ind_to_check.push_back(j);
    }
  }

  return ind_to_check;
}

// Clumping within a distance in bp
// [[Rcpp::export]]
LogicalVector clumping_chr(Environment BM,
                           const IntegerVector& rowInd,
                           const IntegerVector& colInd,
                           const IntegerVector& ordInd,
                           const NumericVector& pos,
                           const NumericVector& sumX,
                           const NumericVector& denoX,
                           double size,
                           double thr) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  size_t n = macc.nrow();
  size_t m = macc.ncol();

  LogicalVector keep(m, false);
  std::vector<size_t> ind_to_check; ind_to_check.reserve(m);

  for (size_t k = 0; k < m; k++) {

    size_t j0 = ordInd[k] - 1;
    keep[j0] = true;

    ind_to_check = which_to_check(j0, keep, pos, size, ind_to_check);

    for (auto&& j : ind_to_check) {

      double xySum = 0;
      for (size_t i = 0; i < n; i++) {
        xySum += macc(i, j) * macc(i, j0);
      }
      double num = xySum - sumX[j] * sumX[j0] / n;
      double r2 = num * num / (denoX[j] * denoX[j0]);
      if (r2 > thr) {
        keep[j0] = false;  // prune
        break;
      }
    }
  }

  return keep;
}

/******************************************************************************/
