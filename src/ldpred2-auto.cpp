/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/utils.h>
#include <bigsparser/SFBM.h>

/******************************************************************************/

inline double square(double x) {
  return x * x;
}

/******************************************************************************/

// [[Rcpp::export]]
List ldpred2_gibbs_auto(Environment corr,
                        const NumericVector& beta_hat,
                        const NumericVector& beta_init,
                        const IntegerVector& order,
                        const NumericVector& n_vec,
                        double p_init,
                        double h2_init,
                        int burn_in,
                        int num_iter,
                        bool verbose = false) {

  XPtr<SFBM> sfbm = corr["address"];

  int m = beta_hat.size();
  myassert_size(sfbm->nrow(), m);
  myassert_size(sfbm->ncol(), m);
  myassert_size(order.size(), m);
  myassert_size(beta_init.size(), m);
  myassert_size(n_vec.size(), m);

  std::vector<double> is_causal(m, p_init);
  arma::vec curr_beta(beta_init.begin(), m);
  arma::vec avg_beta(m, arma::fill::zeros), avg_postp(m, arma::fill::zeros);

  int num_iter_tot = burn_in + num_iter;
  std::vector<double> p_est(num_iter_tot), h2_est(num_iter_tot);

  double nb_causal = m * p_init;
  double cur_h2_est = arma::dot(curr_beta, sfbm->prod(curr_beta));
  double p = p_init, h2 = h2_init, avg_p = 0, avg_h2 = 0;

  for (int k = 0; k < num_iter_tot; k++) {

    for (const int& j : order) {

      double dotprod = sfbm->dot_col(j, curr_beta);
      double res_beta_hat_j = beta_hat[j] + curr_beta[j] - dotprod;

      double C1 = h2 * n_vec[j] / (m * p);
      double C2 = 1 / (1 + 1 / C1);
      double C3 = C2 * res_beta_hat_j ;
      double C4 = ::sqrt(C2 / n_vec[j]);

      double postp = 1 /
        (1 + (1 - p) / p * ::sqrt(1 + C1) * ::exp(-square(C3 / C4) / 2));

      if (k >= burn_in) {
        avg_postp[j] += postp;
        avg_beta[j]  += C3 * postp;
      }

      double nb_rm = is_causal[j];
      is_causal[j] = postp > ::unif_rand();
      nb_causal += is_causal[j] - nb_rm;

      double prev_beta = curr_beta[j];
      curr_beta[j] = is_causal[j] ? (C3 + ::norm_rand() * C4) : 0;
      double diff = curr_beta[j] - prev_beta;
      cur_h2_est += diff * (2 * dotprod + diff);
    }

    p = ::Rf_rbeta(1 + nb_causal, 1 + m - nb_causal);
    h2 = std::max(cur_h2_est, 1e-4);
    if (verbose) Rcout << k + 1 << ": " << p << " // " << h2 << std::endl;

    if (k >= burn_in) {
      avg_p  += p;
      avg_h2 += h2;
    }
    p_est[k]  = p;
    h2_est[k] = h2;
  }

  double est_p  = avg_p  / num_iter;
  double est_h2 = avg_h2 / num_iter;
  if (verbose) Rcout << "Overall: " << est_p << " // " << est_h2 << std::endl;

  return List::create(
    _["beta_est"]    = avg_beta  / num_iter,
    _["postp_est"]   = avg_postp / num_iter,
    _["p_est"]       = est_p,
    _["h2_est"]      = est_h2,
    _["path_p_est"]  = p_est,
    _["path_h2_est"] = h2_est);
}

/******************************************************************************/
