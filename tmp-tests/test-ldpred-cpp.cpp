// [[Rcpp::plugins(cpp11)]]
#define ARMA_64BIT_WORD

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix ldpred_gibbs(const arma::sp_mat& corr,
                           const NumericVector& betas_hat,
                           const NumericVector& n_vec,
                           const NumericVector& p_vec,
                           double coeff,
                           int burn_in = 10,
                           int num_iter = 60,
                           bool sparse = false) {

  int m = betas_hat.size();
  int np = p_vec.size();
  NumericMatrix ldpred_shrink(m, np);

  #pragma omp parallel for num_threads(4)
  for (int i = 0; i < np; i++) {

    double p = p_vec[i];
    #pragma omp critical
    Rcout << "p = " << p << std::endl;

    // some useful constants
    double L = coeff / p;
    double C1 = (1 - p) / p * sqrt(1 + L);
    double C2 = L / (L + 1);
    NumericVector C3 = (-C2 / 2) * n_vec;
    NumericVector C4 = sqrt(C2 / n_vec);

    NumericVector curr_post_means(m), avg_betas(m);
    arma::vec curr_betas(m, arma::fill::zeros);
    // for (int j = 0; j < m; j++) curr_betas[j] = betas_hat[j];

    for (int k = 1; k <= num_iter; k++) {

      // print(k)
      // print(h2_est <- max(0.00001, crossprod(curr_betas)))
      // alpha <- 1 #min(0.99, 1 / h2_est, (h2 + 1 / sqrt(N)) / h2_est)

      for (int j = 0; j < m; j++) {
        curr_betas[j] = 0;
        double res_beta_hat_j = betas_hat[j] - arma::dot(corr.col(j), curr_betas);
        double postp = 1 / (1 + C1 * ::exp(C3[j] * res_beta_hat_j * res_beta_hat_j));
        if (sparse && (postp < p)) {
          curr_betas[j] = curr_post_means[j] = 0;
        } else {
          curr_post_means[j] = C2 * postp * res_beta_hat_j;
          curr_betas[j] = (postp > ::unif_rand())
            ? C4[j] * ::norm_rand() + C2 * res_beta_hat_j : 0;
        }
      }

      if (k > burn_in) avg_betas += curr_post_means;
    }

    #pragma omp critical
    ldpred_shrink(_, i) = (avg_betas / (num_iter - burn_in)) / betas_hat;
  }

  return ldpred_shrink;
}
