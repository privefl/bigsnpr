/******************************************************************************/

// [[Rcpp::plugins(cpp11)]]
#define ARMA_64BIT_WORD

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/utils.h>

using namespace Rcpp;

/******************************************************************************/

inline double square(double x) {
  return x * x;
}

/******************************************************************************/

List ldpred2_gibbs_auto_one(const arma::sp_mat& corr,
                            const NumericVector& betas_hat,
                            const NumericVector& n_vec,
                            double h2_init,
                            double p_init,
                            int burn_in,
                            int num_iter,
                            double h2_min = 1e-4,
                            double h2_max = 1,
                            double prob_jump_to_0 = 1e-4) {

  int m = betas_hat.size();
  std::vector<double> is_causal(m, p_init);
  arma::vec curr_betas(m); curr_betas.fill(::sqrt(h2_init / m));
  double cur_h2_est = arma::dot(curr_betas, corr * curr_betas);

  arma::vec curr_post_means(m, arma::fill::zeros);
  arma::vec avg_betas(m, arma::fill::zeros);
  IntegerVector random_order(m);

  int num_iter_tot = burn_in + num_iter;
  std::vector<double> p_est(num_iter_tot), h2_est(num_iter_tot);

  double p = p_init, h2 = h2_init, alpha = 1;
  double nb_causal = m * p_init;
  double avg_p = 0, avg_h2 = 0;
  int c = 0;

  for (int k = 0; k < num_iter_tot; k++) {

    #pragma omp critical
    random_order = sample(m, m, false, R_NilValue, false);

    for (const int& j : random_order) {

      double dotprod = arma::dot(corr.col(j), curr_betas);
      double res_beta_hat_j = betas_hat[j] + curr_betas[j] - dotprod;

      double C1 = h2 * n_vec[j] / (m * p);
      double C2 = 1 / (1 + 1 / C1);
      double C3 = C2 * res_beta_hat_j ;
      double C4 = ::sqrt(C2 / n_vec[j]);

      double nb_rm = is_causal[j];
      double postp = alpha /
        (1 + (1 - p) / p * ::sqrt(1 + C1) * ::exp(-0.5 * square(C3 / C4)));

      is_causal[j] = postp > ::unif_rand();
      nb_causal += is_causal[j] - nb_rm;
      double prev_beta = curr_betas[j];
      curr_betas[j] = is_causal[j] ? (C3 + ::norm_rand() * C4) : 0;
      double diff = curr_betas[j] - prev_beta;
      cur_h2_est += diff * (2 * dotprod + diff);
      curr_post_means[j] = C3 * postp;
    }

    p = ::Rf_rbeta(1 + nb_causal, 1 + m - nb_causal);
    h2 = std::max(cur_h2_est, h2_min);
    alpha = std::min(1 - prob_jump_to_0, h2_max / h2);

    if (k >= burn_in) {
      c++;
      avg_betas += curr_post_means;
      avg_p     += p;
      avg_h2    += h2;
    }
    Rcout << k + 1 << ": " << p << " // " << h2 << std::endl;
    p_est[k] = p;
    h2_est[k] = h2;
  }

  Rcout << (avg_p / c) << " // " << (avg_h2 / c) << std::endl;

  return List::create(
    _["beta_est"]   = avg_betas / c,
    _["p_est"]      = avg_p     / c,
    _["h2_est"]     = avg_h2    / c,
    _["vec_p_est"]  = wrap(p_est),
    _["vec_h2_est"] = wrap(h2_est));
}

/******************************************************************************/

// [[Rcpp::export]]
List ldpred2_gibbs_auto(const arma::sp_mat& corr,
                        const NumericVector& betas_hat,
                        const NumericVector& n_vec,
                        const NumericVector& h2_init,
                        const NumericVector& p_init,
                        int burn_in,
                        int num_iter,
                        double h2_min = 1e-4,
                        double h2_max = 1,
                        double prob_jump_to_0 = 1e-4,
                        int ncores = 1) {

  int m = betas_hat.size();
  myassert_size(corr.n_rows,  m);
  myassert_size(corr.n_cols,  m);
  myassert_size(n_vec.size(), m);

  int K = p_init.size();
  myassert_size(h2_init.size(), K);

  List res(K);

  #pragma omp parallel for schedule(dynamic, 1) num_threads(ncores)
  for (int k = 0; k < K; k++) {
    List res_k = ldpred2_gibbs_auto_one(
      corr, betas_hat, n_vec, h2_init[k], p_init[k], burn_in, num_iter,
      h2_min, h2_max, prob_jump_to_0);
    #pragma omp critical
    res[k] = res_k;
  }

  return res;
}

/******************************************************************************/
