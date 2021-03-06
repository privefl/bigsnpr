/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/utils.h>

/******************************************************************************/

inline double square(double x) {
  return x * x;
}

/******************************************************************************/

List ldpred2_gibbs_auto_one(const arma::sp_mat& corr,
                            const NumericVector& beta_hat,
                            const NumericVector& beta_init,
                            const NumericVector& order,
                            const NumericVector& n_vec,
                            double p_init,
                            int burn_in,
                            int num_iter,
                            double h2_min,
                            double h2_max,
                            double prob_jump_to_0,
                            bool verbose) {

  int m = beta_hat.size();
  std::vector<double> is_causal(m, p_init);
  arma::vec curr_beta(beta_init.begin(), m);
  arma::vec post_mean_beta(m);
  arma::vec avg_beta(m, arma::fill::zeros);

  int num_iter_tot = burn_in + num_iter;
  std::vector<double> p_est(num_iter_tot), h2_est(num_iter_tot);

  double nb_causal = m * p_init;
  double cur_h2_est = arma::dot(curr_beta, corr * curr_beta);
  double p = p_init, h2 = cur_h2_est, alpha = 1, avg_p = 0, avg_h2 = 0;

  for (int k = 0; k < num_iter_tot; k++) {

    for (const int& j : order) {

      double dotprod = arma::dot(corr.col(j), curr_beta);
      double res_beta_hat_j = beta_hat[j] + curr_beta[j] - dotprod;

      double C1 = h2 * n_vec[j] / (m * p);
      double C2 = 1 / (1 + 1 / C1);
      double C3 = C2 * res_beta_hat_j ;
      double C4 = ::sqrt(C2 / n_vec[j]);

      double postp = alpha /
        (1 + (1 - p) / p * ::sqrt(1 + C1) * ::exp(-square(C3 / C4) / 2));
      post_mean_beta[j] = C3 * postp;

      double nb_rm = is_causal[j];
      is_causal[j] = postp > ::unif_rand();
      nb_causal += is_causal[j] - nb_rm;

      double prev_beta = curr_beta[j];
      curr_beta[j] = is_causal[j] ? (C3 + ::norm_rand() * C4) : 0;
      double diff = curr_beta[j] - prev_beta;
      cur_h2_est += diff * (2 * dotprod + diff);
    }

    p = ::Rf_rbeta(1 + nb_causal, 1 + m - nb_causal);
    h2 = std::max(cur_h2_est, h2_min);
    alpha = std::min(1 - prob_jump_to_0, h2_max / h2);

    if (k >= burn_in) {
      avg_beta += post_mean_beta;
      avg_p    += p;
      avg_h2   += h2;
    }
    p_est[k] = p;
    h2_est[k] = h2;
    if (verbose) Rcout << k + 1 << ": " << p << " // " << h2 << std::endl;
  }

  double est_p  = avg_p  / num_iter;
  double est_h2 = avg_h2 / num_iter;
  if (verbose) Rcout << "Overall: " << est_p << " // " << est_h2 << std::endl;

  return List::create(
    _["beta_est"]    = avg_beta / num_iter,
    _["p_est"]       = est_p,
    _["h2_est"]      = est_h2,
    _["path_p_est"]  = p_est,
    _["path_h2_est"] = h2_est);
}

/******************************************************************************/

// [[Rcpp::export]]
List ldpred2_gibbs_auto(const arma::sp_mat& corr,
                        const NumericVector& beta_hat,
                        const NumericVector& beta_init,
                        const NumericVector& order,
                        const NumericVector& n_vec,
                        const NumericVector& p_init,
                        int burn_in,
                        int num_iter,
                        double h2_min = 1e-4,
                        double h2_max = 1,
                        double prob_jump_to_0 = 1e-4,
                        bool verbose = false,
                        int ncores = 1) {

  int m = beta_hat.size();
  myassert_size(corr.n_rows, m);
  myassert_size(corr.n_cols, m);
  myassert_size(order.size(), m);
  myassert_size(beta_init.size(), m);
  myassert_size(n_vec.size(), m);

  int K = p_init.size();
  List res(K);

  #pragma omp parallel for schedule(dynamic, 1) num_threads(ncores)
  for (int k = 0; k < K; k++) {

    List res_k = ldpred2_gibbs_auto_one(
      corr, beta_hat, beta_init, order, n_vec, p_init[k], burn_in, num_iter,
      h2_min, h2_max, prob_jump_to_0, verbose);

    #pragma omp critical
    res[k] = res_k;
  }

  return res;
}

/******************************************************************************/
