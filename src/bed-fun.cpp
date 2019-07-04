/******************************************************************************/

#include "bed-acc.h"
#include <bigstatsr/colstats.hpp>

/******************************************************************************/

// [[Rcpp::export]]
ListOf<NumericVector> bedcolvars(const std::string path,
                                 int n_total, int m_total,
                                 const IntegerVector& row_ind,
                                 const IntegerVector& col_ind,
                                 const RawMatrix& lookup_byte) {

  bedAcc macc(path, n_total, m_total, row_ind, col_ind, lookup_byte);
  size_t n = macc.nrow();
  size_t m = macc.ncol();

  NumericVector res(m), res2(m), res3(m);
  double x, xSum, xxSum, c;

  for (size_t j = 0; j < m; j++) {
    xSum = xxSum = c = 0;
    c = 0;
    for (size_t i = 0; i < n; i++) {
      x = macc(i, j);
      if (x != 3) {
        xSum += x;
        xxSum += x*x;
        c++;
      }
    }
    res[j] = xxSum - xSum * xSum / n;
    res2[j] = xSum;
    res3[j] = c;
  }

  return List::create(_["sum"]  = res2,
                      _["var"]  = res / (res3 - 1),
                      _["nona"] = res3);
}

/******************************************************************************/

// Clumping within a distance in bp (directly on a bed file)
// [[Rcpp::export]]
LogicalVector bed_clumping_chr(const std::string path, int n_total, int m_total,
                               const IntegerVector& row_ind,
                               const IntegerVector& col_ind,
                               const RawMatrix& lookup_byte,
                               const NumericMatrix& lookup_scale,
                               const IntegerVector& ordInd,
                               const IntegerVector& pos,
                               int size,
                               double thr) {

  bedAccScaled macc(path, n_total, m_total, row_ind, col_ind, lookup_byte,
                    lookup_scale);

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
