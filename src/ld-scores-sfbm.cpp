/******************************************************************************/

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include <bigsparser/SFBM.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
NumericVector ld_scores_sfbm(Rcpp::Environment X, bool compact) {

  Rcpp::XPtr<SFBM> sfbm = X["address"];
  const NumericVector p = X["p"];
  int m = sfbm->ncol();
  const double * data = sfbm->i_x();

  NumericVector res(m);

  if (compact) {

    for (int j = 0; j < m; j++) {

      size_t lo = p[j];
      size_t up = p[j + 1];

      for (size_t k = lo; k != up; k++)
        res[j] += data[k] * data[k];
    }

  } else {

    for (int j = 0; j < m; j++) {

      size_t lo = 2 * p[j];
      size_t up = 2 * p[j + 1];

      for (size_t k = lo + 1; k < up; k += 2)
        res[j] += data[k] * data[k];
    }

  }

  return res;
}

/******************************************************************************/
