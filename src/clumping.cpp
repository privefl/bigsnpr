/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <RcppParallel.h>

using namespace RcppParallel;

/******************************************************************************/

std::vector<int>& which_to_check(int j0,
                                 const LogicalVector& keep,
                                 const NumericVector& pos,
                                 double size,
                                 std::vector<int>& ind_to_check) {

  ind_to_check.clear();
  int m = pos.size();

  int pos_min = pos[j0] - size;
  int pos_max = pos[j0] + size;
  bool not_min = true;
  bool not_max = true;

  for (int l = 1; not_max || not_min; l++) {
    if (not_max) {
      int j = j0 + l;
      not_max = (j < m) && (pos[j] <= pos_max);   // within a window..
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

struct Prune : public Worker {

  SubBMCode256Acc macc;
  size_t j0;
  std::vector<int> ind_to_check;
  RVector<int> keep;
  RVector<double> sumX, denoX;
  double thr;

  // constructors
  Prune(SubBMCode256Acc macc, LogicalVector& keep,
        const NumericVector& sumX, const NumericVector& denoX, double thr) :
    macc(macc), keep(keep), sumX(sumX), denoX(denoX), thr(thr) {}

  Prune(const Prune& prune, size_t j0, std::vector<int>& ind_to_check) :
    macc(prune.macc), j0(j0), ind_to_check(ind_to_check), keep(prune.keep),
    sumX(prune.sumX), denoX(prune.denoX), thr(prune.thr) {}

  // parallel computations
  void operator()(size_t begin, size_t end) {

    size_t n = macc.nrow();
    for (size_t k = begin; k < end; k++) {

      if (!keep[j0]) break;
      size_t j = ind_to_check[k];

      double xySum = 0;
      for (size_t i = 0; i < n; i++) {
        xySum += macc(i, j) * macc(i, j0);
      }
      double num = xySum - sumX[j] * sumX[j0] / n;
      double r2 = num * num / (denoX[j] * denoX[j0]);

      if (r2 > thr) keep[j0] = false;  // prune
    }
  }
};

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

  size_t m = macc.ncol();

  LogicalVector keep(m, false);
  std::vector<int> ind_to_check; ind_to_check.reserve(m);

  Prune prune_init(macc, keep, sumX, denoX, thr);

  for (size_t k = 0; k < m; k++) {

    size_t j0 = ordInd[k] - 1;
    keep[j0] = true;

    ind_to_check = which_to_check(j0, keep, pos, size, ind_to_check);

    Prune prune(prune_init, j0, ind_to_check);
    parallelFor(0, ind_to_check.size(), prune);
  }

  return keep;
}

/******************************************************************************/
