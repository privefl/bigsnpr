// [[Rcpp::plugins(cpp11)]]
#define ARMA_64BIT_WORD

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ldpred_gibbs_auto(const arma::sp_mat& corr,
                                const NumericVector& betas_hat,
                                const NumericVector& n_vec,
                                double coeff,
                                double p_init,
                                const NumericVector& w,
                                int burn_in = 10,
                                int num_iter = 60,
                                bool sparse = false) {

  int m = betas_hat.size();
  NumericVector postp(m, p_init);

  NumericVector curr_post_means(m), avg_betas(m);
  arma::vec curr_betas(m, arma::fill::zeros);

  double p = p_init, avg_p = 0, avg_h2 = 0;

  for (int k = 1; k <= num_iter; k++) {

    // some useful constants
    double L = coeff / p;
    double C1 = (1 - p) / p * sqrt(1 + L);
    double C2 = L / (L + 1);
    NumericVector C3 = (-C2 / 2) * n_vec;
    NumericVector C4 = sqrt(C2 / n_vec);

    for (int j = 0; j < m; j++) {
      curr_betas[j] = 0;
      double res_beta_hat_j = betas_hat[j] - arma::dot(corr.col(j), curr_betas);
      postp[j] = 1 / (1 + C1 * ::exp(C3[j] * res_beta_hat_j * res_beta_hat_j));
      if (sparse && (postp[j] < p)) {
        curr_betas[j] = curr_post_means[j] = 0;
      } else {
        curr_post_means[j] = C2 * postp[j] * res_beta_hat_j;
        curr_betas[j] = (postp[j] > ::unif_rand())
          ? C4[j] * ::norm_rand() + C2 * res_beta_hat_j : 0;
      }
    }

    p = sum(postp * w) / sum(w);
    double h2 = arma::dot(curr_betas, curr_betas);

    Rcout << p << " // " << h2 << std::endl;

    if (k > burn_in) {
      avg_betas += curr_post_means;
      avg_p += p;
      avg_h2 += h2;
    }
  }

  Rcout << avg_p / (num_iter - burn_in) << " // " <<
    avg_h2 / (num_iter - burn_in) << std::endl;

  return (avg_betas / (num_iter - burn_in)) / betas_hat;
}
