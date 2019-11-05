/******************************************************************************/

#include "bed-acc.h"

/******************************************************************************/

// Clumping within a distance in bp (directly on a bed file)
// [[Rcpp::export]]
LogicalVector bed_clumping_chr(Environment obj_bed,
                               const IntegerVector& ind_row,
                               const IntegerVector& ind_col,
                               const NumericVector& center,
                               const NumericVector& scale,
                               const IntegerVector& ordInd,
                               const IntegerVector& pos,
                               int size,
                               double thr) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAccScaled macc(xp_bed, ind_row, ind_col, center, scale);

  int n = macc.nrow();
  int m = macc.ncol();

  LogicalVector keep(m, false);

  for (int k = 0; k < m; k++) {

    int j0 = ordInd[k] - 1;
    keep[j0] = true;
    bool not_min = true;
    bool not_max = true;
    int pos_min = pos[j0] - size;
    int pos_max = pos[j0] + size;

    for (int l = 1; not_max || not_min; l++) {

      if (not_max) {
        int j = j0 + l;
        not_max = (j < m) && (pos[j] <= pos_max);  // within a window..
        if (not_max && keep[j]) {  // look only at already selected ones
          double r = 0;
          for (int i = 0; i < n; i++) {
            r += macc(i, j) * macc(i, j0);
          }
          if ((r * r) > thr) {
            // Rcout << r << std::endl;
            keep[j0] = false;  // prune
            break;
          }
        }
      }

      if (not_min) {
        int j = j0 - l;
        not_min = (j >= 0) && (pos[j] >= pos_min);  // within a window..
        if (not_min && keep[j]) {  // look only at already selected ones
          double r = 0;
          for (int i = 0; i < n; i++) {
            r += macc(i, j) * macc(i, j0);
          }
          if ((r * r) > thr) {
            keep[j0] = false;  // prune
            break;
          }
        }
      }
    }
  }

  return keep;
}

/******************************************************************************/
