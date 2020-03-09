/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include "clumping-utils.h"

/******************************************************************************/

// Clumping within a distance
// [[Rcpp::export]]
void clumping_chr(Environment BM,
                  Environment BM2,
                  const IntegerVector& rowInd,
                  const IntegerVector& colInd,
                  const IntegerVector& ordInd,
                  const IntegerVector& rankInd,
                  const NumericVector& pos,
                  const NumericVector& sumX,
                  const NumericVector& denoX,
                  double size,
                  double thr,
                  int ncores) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  XPtr<FBM_RW> xpBM2 = BM2["address_rw"];
  int * keep = static_cast<int *>(xpBM2->matrix());

  size_t n = macc.nrow();
  size_t m = macc.ncol();
  // Rcout << ncores << std::endl;

  #pragma omp parallel num_threads(ncores)
  {
    std::vector<int> ind_to_check; ind_to_check.reserve(m);

    #pragma omp for schedule(dynamic, 1)
    for (size_t k = 0; k < m; k++) {

      // Rcout << k << std::endl;
      size_t j0 = ordInd[k] - 1;
      // Rcout << rankInd[j0] << " == " << k + 1 << std::endl;

      // which to check
      ind_to_check = which_to_check(j0, keep, rankInd, pos, size, ind_to_check);
      int nb_check = ind_to_check.size();
      // Rcout << nb_check << std::endl;

      bool keep_j0 = true;
      bool all_checked = false;
      while (!all_checked) {

        all_checked = true;
        for (int k2 = 0; k2 < nb_check; k2++) {

          int jk = ind_to_check[k2];

          // check if already done with this index
          if (jk != -1) {
            if (keep[jk] == -1) {
              // waiting for the thread checking 'j' to finish
              all_checked = false;
            } else if (keep[jk] == 0) {
              // no need to check this one after all
              ind_to_check[k2] = -1;
            } else {
              size_t j = jk;
              double xySum = 0;
              for (size_t i = 0; i < n; i++) {
                xySum += macc(i, j) * macc(i, j0);
              }
              double num = xySum - sumX[j] * sumX[j0] / n;
              double r2 = num * num / (denoX[j] * denoX[j0]);

              if (r2 > thr) {
                keep_j0 = false;  // prune
                all_checked = true;
                break;
              } else {
                ind_to_check[k2] = -1;
              }
            }
          }
        }
      }

      #pragma omp atomic write
      keep[j0] = keep_j0;  // can now be used by other threads
    }
  }
}

/******************************************************************************/
