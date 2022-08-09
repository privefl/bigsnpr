/******************************************************************************/

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include <bigsparser/SFBM.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
NumericVector ld_scores_sfbm(Rcpp::Environment X, bool compact, int ncores) {

  Rcpp::XPtr<SFBM> sfbm = X["address"];
  const NumericVector p = X["p"];
  int m = sfbm->ncol();
  const double * data = sfbm->i_x();

  std::vector<double> res(m);

  int chunk_size = ceil(m / (10.0 * ncores));

  if (compact) {

    #pragma omp parallel for schedule(dynamic, chunk_size) num_threads(ncores)
    for (int j = 0; j < m; j++) {

      size_t lo = p[j];
      size_t up = p[j + 1];

      double ld_j = 0;

      for (size_t k = lo; k != up; k++)
        ld_j += data[k] * data[k];

      res[j] = ld_j;
    }

  } else {

    #pragma omp parallel for schedule(dynamic, chunk_size) num_threads(ncores)
    for (int j = 0; j < m; j++) {

      size_t lo = 2 * p[j];
      size_t up = 2 * p[j + 1];

      double ld_j = 0;

      for (size_t k = lo + 1; k < up; k += 2)
        ld_j += data[k] * data[k];

      res[j] = ld_j;
    }

  }

  return Rcpp::wrap(res);
}

/******************************************************************************/
