#ifndef CLUMPING_UTILS_HPP
#define CLUMPING_UTILS_HPP

/******************************************************************************/

#include <Rcpp.h>

using namespace Rcpp;

/******************************************************************************/

inline std::vector<int>& which_to_check(int j0,
                                        const int * keep,
                                        const IntegerVector& rankInd,
                                        const NumericVector& pos,
                                        double size,
                                        std::vector<int>& ind_to_check) {

  ind_to_check.clear();
  int m = pos.size();

  int pos_min = pos[j0] - size;
  int pos_max = pos[j0] + size;
  bool not_min = true;
  bool not_max = true;

  for (int k = 1; not_max || not_min; k++) {
    if (not_max) {
      int j = j0 + k;
      not_max = (j < m) && (pos[j] <= pos_max);   // within a window..
      if (not_max && (rankInd[j0] > rankInd[j]) && (keep[j] != 0))
        ind_to_check.push_back(j);
    }
    if (not_min) {
      int j = j0 - k;
      not_min = (j >= 0) && (pos[j] >= pos_min);  // within a window..
      if (not_min && (rankInd[j0] > rankInd[j]) && (keep[j] != 0))
        ind_to_check.push_back(j);
    }
  }

  return ind_to_check;
}

/******************************************************************************/

#endif // CLUMPING_UTILS_HPP
