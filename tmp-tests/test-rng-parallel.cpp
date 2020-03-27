#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double test_rng_openmp(int nb, int ncores) {

  double sum = 0;

  #pragma omp parallel for reduction(+:sum) num_threads(ncores)
  for (int k = 0; k < nb; k++) {
    sum += ::sqrt(::unif_rand());
    sum += ::exp(::norm_rand());
  }

  return sum;
}

// [[Rcpp::export]]
double test_rng_openmp2(int nb, int m, int ncores) {

  double sum = 0;
  IntegerVector random_order(m);

  #pragma omp parallel for reduction(+:sum) num_threads(ncores)
  for (int k = 0; k < nb; k++) {
    #pragma omp critical
    random_order = sample(m, m, false, R_NilValue, false);
    sum += random_order[0];
  }

  return sum;
}

/*** R
# test_rng_openmp(2e9, 4)
test_rng_openmp2(1e7, 100, 4)
*/
