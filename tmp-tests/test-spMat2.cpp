#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::sp_mat test_arma(int n = 5000, double p = 0.01) {

  arma::sp_mat mat(n, n);

  #pragma omp parallel for num_threads(4)
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      if (unif_rand() < p) {
        #pragma omp critical
        mat(i, j) = unif_rand();
      }

  return mat;
}

/*** R
library(Matrix)
sp <- forceSymmetric(test_arma())
sp_as_list <- lapply(1:ncol(sp), function(i) {
  print(i)
  ind <- which(sp[, i] != 0)
  list(ind, sp[ind, i])
})
*/
