// [[Rcpp::plugins(cpp11)]]
#define ARMA_64BIT_WORD

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ldpred_gibbs_auto3(const arma::sp_mat& corr,
                                 const NumericVector& betas_hat,
                                 const IntegerVector& order,
                                 const NumericVector& n_vec,
                                 double h2_init,
                                 double p_init,
                                 const NumericVector& w,
                                 int burn_in = 10,
                                 int num_iter = 60,
                                 bool sparse = false) {

  int m = betas_hat.size();
  NumericVector post_p(m, p_init);
  NumericVector prop_w = w / sum(w);

  NumericVector curr_post_means(m), avg_betas(m);
  arma::vec curr_betas(m, arma::fill::zeros);
  // std::copy(betas_hat.begin(), betas_hat.end(), curr_betas.begin());

  int num_iter_tot = burn_in + num_iter;
  NumericVector p_est(num_iter_tot), h2_est(num_iter_tot);

  double p = p_init, h2 = h2_init, avg_p = 0, avg_h2 = 0;

  for (int k = 0; k < num_iter_tot; k++) {

    // double p = sum(post_p * prop_w);
    // if (p < 1e-5) p = 1e-5;

    for (const int& j : order) {

      curr_betas[j] = 0;
      double res_beta_hat_j = betas_hat[j] - arma::dot(corr.col(j), curr_betas);

      double L = h2 * n_vec[j] / (m * p);
      // if (j == 1) Rcout << L << std::endl;
      double C2 = 1 / (1 + 1 / L);
      double prev_post_p = post_p[j];
      post_p[j] = 1 / (1 + (1 - p) / p * ::sqrt(1 + L) *
        ::exp(-C2 * n_vec[j] / 2 * res_beta_hat_j * res_beta_hat_j));
      p += (post_p[j] - prev_post_p) * prop_w[j];

      if (sparse && (post_p[j] < p)) {
        curr_betas[j] = curr_post_means[j] = 0;
      } else {
        curr_post_means[j] = C2 * post_p[j] * res_beta_hat_j;
        curr_betas[j] = (post_p[j] > ::unif_rand())
          ? ::sqrt(C2 / n_vec[j]) * ::norm_rand() + C2 * res_beta_hat_j : 0;
      }

    }

    h2 = arma::dot(curr_betas, curr_betas);
    if (k >= burn_in) {
      avg_betas += curr_post_means;
      avg_p += p;
      avg_h2 += h2;
    }
    Rcout << p << " // " << h2 << std::endl;
    p_est[k] = p;
    h2_est[k] = h2;
  }

  Rcout << (avg_p / num_iter) << " // " << (avg_h2 / num_iter) << std::endl;

  return List::create(
    _["shrink"] = (avg_betas / num_iter) / betas_hat,
    _["p_est"]  = avg_p / num_iter,
    _["h2_est"] = avg_h2 / num_iter,
    _["vec_p_est"]  = p_est,
    _["vec_h2_est"] = h2_est);
}
