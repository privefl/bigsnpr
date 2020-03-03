/******************************************************************************/

#include "bed-acc.h"
#include "clumping-utils.h"

/******************************************************************************/

struct Prune : public Worker {

  bedAccScaled macc;
  size_t j0;
  std::vector<int> ind_to_check;
  RVector<int> keep;
  double thr;

  // constructors
  Prune(bedAccScaled macc, LogicalVector& keep, double thr) :
    macc(macc), keep(keep), thr(thr) {}

  Prune(const Prune& prune, size_t j0, std::vector<int>& ind_to_check) :
    macc(prune.macc), j0(j0), ind_to_check(ind_to_check), keep(prune.keep),
    thr(prune.thr) {}

  // parallel computations
  void operator()(size_t begin, size_t end) {

    size_t n = macc.nrow();
    for (size_t k = begin; k < end; k++) {

      if (!keep[j0]) break;
      size_t j = ind_to_check[k];

      double r = 0;
      for (size_t i = 0; i < n; i++) {
        r += macc(i, j) * macc(i, j0);
      }
      double r2 = r * r;

      if (r2 > thr) keep[j0] = false;  // prune
    }
  }
};

/******************************************************************************/

// Clumping within a distance (directly on a bed file)
// [[Rcpp::export]]
LogicalVector bed_clumping_chr(Environment obj_bed,
                               const IntegerVector& ind_row,
                               const IntegerVector& ind_col,
                               const NumericVector& center,
                               const NumericVector& scale,
                               const IntegerVector& ordInd,
                               const NumericVector& pos,
                               double size,
                               double thr) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAccScaled macc(xp_bed, ind_row, ind_col, center, scale);

  int m = macc.ncol();
  LogicalVector keep(m, false);
  std::vector<int> ind_to_check; ind_to_check.reserve(m);

  Prune prune_init(macc, keep, thr);

  for (int k = 0; k < m; k++) {

    size_t j0 = ordInd[k] - 1;
    keep[j0] = true;

    ind_to_check = which_to_check(j0, keep, pos, size, ind_to_check);

    Prune prune(prune_init, j0, ind_to_check);
    parallelFor(0, ind_to_check.size(), prune);
  }

  return keep;
}

/******************************************************************************/
