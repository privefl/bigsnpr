/******************************************************************************/

// [[Rcpp::depends(RcppParallel, rmio, bigstatsr)]]
#include <bigstatsr/BMCodeAcc.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

/******************************************************************************/

int which_pos_min(const NumericVector& pos, int j0, double size) {

  double pos_min = pos[j0] - size;
  int begin = 0, end = j0;

  while (true) {

    int mid = (begin + end) / 2;

    if (mid == begin) {
      if (pos[begin] < pos_min) {
        return end;
      } else {
        return begin;
      }
    } else {
      if (pos[mid] < pos_min) {
        begin = mid;
      } else {
        end = mid;
      }
    }
  }
}

int which_pos_max(const NumericVector& pos, int j0, double size) {

  double pos_max = pos[j0] + size;
  int begin = j0, end = pos.size() - 1;

  while (true) {

    int mid = (begin + end) / 2;

    if (mid == begin) {
      if (pos[end] > pos_max) {
        return begin;
      } else {
        return end;
      }
    } else {
      if (pos[mid] > pos_max) {
        end = mid;
      } else {
        begin = mid;
      }
    }
  }
}

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

    for (size_t j = begin; j < end; j++) {
      if (remain[j]) { // if already excluded, goto next
        double xySum = 0;
        for (size_t i = 0; i < n; i++) {
          xySum += macc(i, j) * macc(i, j0);
        }
        double num = xySum - sumX[j] * sumX[j0] / n;
        double r2 = num * num / (denoX[j] * denoX[j0]);
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
                       const NumericVector& pos,
                       LogicalVector& remain,
                       const NumericVector& sumX,
                       const NumericVector& denoX,
                       double size,
                       double thr) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);
  int m = macc.ncol();

  LogicalVector keep(m); // init with all false

  Prune prune_init(macc, remain, sumX, denoX, thr);

  int grain = 1; //::sqrt(size);

  for (int k = 0; k < m; k++) {
    int j0 = ordInd[k] - 1;
    if (remain[j0]) { // if already excluded, goto next
      remain[j0] = false;
      keep[j0] = true;
      Prune prune(prune_init, j0);
      int j_min = which_pos_min(pos, j0, size);
      int j_max = which_pos_max(pos, j0, size);
      parallelFor(j_min, j_max, prune, grain);
    }
  }

  return keep;
}

/******************************************************************************/
