/******************************************************************************/

#include <bigsparser/SFBM.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
NumericVector ld_scores_sfbm(Rcpp::Environment X,
                             const IntegerVector& ind_sub,
                             int ncores) {

  Rcpp::XPtr<SFBM> sfbm = X["address"];
  const NumericVector p = X["p"];
  const double * data = sfbm->i_x();

  int m2 = sfbm->ncol();
  std::vector<bool> use(m2, false);
  for (const int& k : ind_sub) use[k] = true;

  int m = ind_sub.size();
  std::vector<double> res(m);

  int chunk_size = ceil(m / (10.0 * ncores));

  if (sfbm->is_compact()) {

    std::vector<int> first_i = X["first_i"];

    #pragma omp parallel for schedule(dynamic, chunk_size) num_threads(ncores)
    for (int j = 0; j < m; j++) {

      int j2 = ind_sub[j];
      size_t lo = p[j2];
      size_t up = p[j2 + 1];
      int i = first_i[j2];

      double ld_j = 0;

      for (size_t k = lo; k != up; k++, i++)
        if (use[i])
          ld_j += data[k] * data[k];

      res[j] = ld_j;
    }

  } else {

    #pragma omp parallel for schedule(dynamic, chunk_size) num_threads(ncores)
    for (int j = 0; j < m; j++) {

      int j2 = ind_sub[j];
      size_t lo = 2 * p[j2];
      size_t up = 2 * p[j2 + 1];

      double ld_j = 0;

      for (size_t k = lo + 1; k < up; k += 2)
        if (use[data[k - 1]])
          ld_j += data[k] * data[k];

      res[j] = ld_j;
    }

  }

  return Rcpp::wrap(res);
}

/******************************************************************************/
