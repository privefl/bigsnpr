#define ARMA_USE_SUPERLU
// [Rcpp::depends(RcppArmadillo)]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::vec sp_solver(arma::sp_mat K, arma::vec x) {
  arma::superlu_opts opts;
  opts.symmetric  = true;
  arma::vec res;
  arma::spsolve(res, K, x, "superlu", opts);
  return res;
}

/*** R
library(Matrix)
K <- sparseMatrix(i = c(1, 2, 1), j = c(1, 2, 2), x = c(1, 1, 0.5), symmetric = TRUE)
x <- runif(2)
sp_solver(K, x)
*/
