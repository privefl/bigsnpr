// [[Rcpp::plugins(cpp11)]]
#define ARMA_64BIT_WORD

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ldpred_gibbs_auto3(const arma::sp_mat& corr,
                                 const NumericVector& betas_hat,
                                 const IntegerVector& order,
                                 const NumericVector& n_vec,
                                 double h2_max,
                                 double p_init,
                                 int burn_in,
                                 int num_iter) {

  int m = betas_hat.size();
  NumericVector post_p(m, p_init), curr_post_means(m), avg_betas(m);
  arma::vec curr_betas(m, arma::fill::zeros);

  int num_iter_tot = burn_in + num_iter;
  NumericVector p_est(num_iter_tot), h2_est(num_iter_tot);

  double p = p_init, h2 = h2_max, alpha = 1;
  double avg_p = 0, avg_h2 = 0;

  int k = 0;
  for (; k < num_iter_tot; k++) {

    for (const int& j : order) {

      curr_betas[j] = 0;
      double res_beta_hat_j = betas_hat[j] - arma::dot(corr.col(j), curr_betas);

      double L = h2 * n_vec[j] / (m * p);
      double C2 = 1 / (1 + 1 / L);

      post_p[j] = 1 / (1 + (1 - p) / p * ::sqrt(1 + L) *
        ::exp(-C2 * n_vec[j] / 2 * res_beta_hat_j * res_beta_hat_j));

      curr_post_means[j] = C2 * post_p[j] * res_beta_hat_j;
      curr_betas[j] = ((alpha * post_p[j]) > ::unif_rand())
        ? C2 * res_beta_hat_j + ::norm_rand() * ::sqrt(C2 / n_vec[j]) : 0;
    }

    double samp_mean_p = Rcpp::mean(sample(post_p, m, true));
    p = std::max(1e-5, samp_mean_p);
    arma::vec samp_betas = Rcpp::RcppArmadillo::sample(curr_betas, m, true);
    h2 = std::max(1e-4, arma::dot(samp_betas, samp_betas));
    alpha = std::min(1.0, h2_max / h2);

    if (k >= burn_in) {
      avg_betas += curr_post_means;
      avg_p += p;
      avg_h2 += h2;
    }
    Rcout << p << " // " << h2 << std::endl;
    p_est[k] = p;
    h2_est[k] = h2;
  }

  Rcout << (avg_p / (k - burn_in)) << " // " <<
    (avg_h2 / (k - burn_in)) << std::endl;

  return List::create(
    _["beta_est"] = (avg_betas / (k - burn_in)),
    _["p_est"]  = avg_p / (k - burn_in),
    _["h2_est"] = avg_h2 / (k - burn_in),
    _["vec_p_est"]  = p_est,
    _["vec_h2_est"] = h2_est);
}
