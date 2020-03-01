/******************************************************************************/

// [[Rcpp::depends(RcppParallel, rmio, bigstatsr)]]
#include <bigstatsr/BMCodeAcc.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

/******************************************************************************/

struct Prune : public Worker {

  SubBMCode256Acc macc;
  size_t j0, n;
  RVector<int> remain;
  RVector<double> sumX, denoX;
  double thr;

  // constructors
  Prune(SubBMCode256Acc macc, LogicalVector& remain,
        const NumericVector& sumX, const NumericVector& denoX, double thr) :
    macc(macc), j0(0), n(macc.nrow()), remain(remain),
    sumX(sumX), denoX(denoX), thr(thr) {}
  Prune(const Prune& prune, size_t j0) :
    macc(prune.macc), j0(j0), n(prune.n), remain(prune.remain),
    sumX(prune.sumX), denoX(prune.denoX), thr(prune.thr) {}

  void operator()(size_t begin, size_t end) {

    double xySum, num, r2;
    size_t i, j;

    for (j = begin; j < end; j++) {
      if (remain[j]) { // if already excluded, goto next
        xySum = 0;
        for (i = 0; i < n; i++) {
          xySum += macc(i, j) * macc(i, j0);
        }
        num = xySum - sumX[j] * sumX[j0] / n;
        r2 = num * num / (denoX[j] * denoX[j0]);
        if (r2 > thr) remain[j] = false; // prune
      }
    }
  }
};

// Clumping within a distance in number of SNPs
// [[Rcpp::export]]
LogicalVector clumping(Environment BM,
                       const IntegerVector& rowInd,
                       const IntegerVector& colInd,
                       const IntegerVector& ordInd,
                       LogicalVector& remain,
                       const NumericVector& sumX,
                       const NumericVector& denoX,
                       int size,
                       double thr) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd - 1, colInd - 1, BM["code256"]);
  // int n = macc.nrow();
  int m = macc.ncol();

  // double xySum, num, r2;
  int j0, k, j_min, j_max;

  LogicalVector keep(m); // init with all false

  Prune prune_init(macc, remain, sumX, denoX, thr);

  int grain = ::sqrt(size);

  for (k = 0; k < m; k++) {
    j0 = ordInd[k] - 1;
    if (remain[j0]) { // if already excluded, goto next
      remain[j0] = false;
      keep[j0] = true;
      Prune prune(prune_init, j0);
      j_min = std::max(0, j0 - size);
      j_max = std::min(m, j0 + size + 1);
      parallelFor(j_min, j_max, prune, grain);
      // for (j = j_min; j < j_max; j++) {
      //   if (remain[j]) { // if already excluded, goto next
      //     xySum = 0;
      //     for (i = 0; i < n; i++) {
      //       xySum += macc(i, j) * macc(i, j0);
      //     }
      //     num = xySum - sumX[j] * sumX[j0] / n;
      //     r2 = num * num / (denoX[j] * denoX[j0]);
      //     if (r2 > thr) remain[j] = false; // prune
      //   }
      // }
    }
  }

  return keep;
}

/******************************************************************************/
