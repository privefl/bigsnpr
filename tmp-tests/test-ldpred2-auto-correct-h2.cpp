/******************************************************************************/

// [[Rcpp::depends(rmio, bigstatsr, bigsparser)]]
#include <bigsparser/SFBM.h>
#include <bigstatsr/utils.h>

/******************************************************************************/

// const double MIN_P  = 1e-5;
const double MIN_H2 = 1e-3;

/******************************************************************************/

// [[Rcpp::export]]
List ldpred2_gibbs_auto(Environment corr,
                        const NumericVector& beta_hat,
                        const NumericVector& n_vec,
                        const IntegerVector& ind_sub,
                        double p_init,
                        double h2_init,
                        int burn_in,
                        int num_iter,
                        bool no_jump_sign,
                        double shrink_corr,
                        const NumericVector& p_bounds,
                        double mean_ld = 1,
                        bool verbose = false) {

  XPtr<SFBM> sfbm = corr["address"];

  int m = beta_hat.size();
  myassert_size(n_vec.size(), m);
  NumericVector curr_beta(m);  // only for the subset
  int m2 = sfbm->ncol();
  NumericVector dotprods(m2);  // for the full corr

  NumericVector avg_beta(m), avg_postp(m), avg_beta_hat(m);

  int num_iter_tot = burn_in + num_iter;
  NumericVector p_est(num_iter_tot, NA_REAL), h2_est(num_iter_tot, NA_REAL);

  double cur_h2_est = 0;
  double h2 = std::max(h2_init, MIN_H2);
  double p = std::min(std::max(p_bounds[0], p_init), p_bounds[1]);

  double gap0 = 2 *
    std::inner_product(beta_hat.begin(), beta_hat.end(), beta_hat.begin(), 0.0);

  std::vector<int> ind_causal;

  for (int k = 0; k < num_iter_tot; k++) {

    double inv_odd_p = (1 - p) / p;
    double sigma2 = h2 / (m * p);
    double gap = 0;

    ind_causal.clear();

    for (int j = 0; j < m; j++) {

      int j2 = ind_sub[j];
      double dotprod = dotprods[j2];
      double res_beta_hat_j = beta_hat[j] - shrink_corr * (dotprod - curr_beta[j]);

      double C1 = sigma2 * n_vec[j];
      double C2 = 1 / (1 + 1 / C1);
      double C3 = C2 * res_beta_hat_j;
      double C4 = C2 / n_vec[j];

      double postp = 1 /
        (1 + inv_odd_p * ::sqrt(1 + C1) * ::exp(-C3 * C3 / C4 / 2));

      double prev_beta = curr_beta[j];
      double dotprod_shrunk = shrink_corr * dotprod + (1 - shrink_corr) * prev_beta;

      if (k >= burn_in) {
        avg_postp[j]    += postp;
        avg_beta[j]     += C3 * postp;
        avg_beta_hat[j] += dotprod_shrunk;
      }

      double diff = -prev_beta;
      if (postp > ::unif_rand()) {

        double samp_beta = ::Rf_rnorm(C3, ::sqrt(C4));

        if (no_jump_sign && (samp_beta * prev_beta) < 0) {
          curr_beta[j] = 0;
        } else {
          curr_beta[j] = samp_beta;
          diff += samp_beta;
          ind_causal.push_back(j);
          gap += samp_beta * samp_beta;
        }

      } else {
        curr_beta[j] = 0;
      }

      if (diff != 0) {
        cur_h2_est += diff * (2 * dotprod_shrunk + diff);
        dotprods = sfbm->incr_mult_col(j2, dotprods, diff);
      }
    }

    if (gap > gap0) {
      avg_beta.fill(NA_REAL); avg_postp.fill(NA_REAL); avg_beta_hat.fill(NA_REAL);
      break;
    }

    int nb_causal = ind_causal.size();
    p = ::Rf_rbeta(1 + nb_causal / mean_ld, 1 + (m - nb_causal) / mean_ld);
    p = std::min(std::max(p_bounds[0], p), p_bounds[1]);
    h2 = std::max(cur_h2_est, MIN_H2);

    if (verbose) Rcout <<
      k + 1 << ": " << p << " // " << h2 << std::endl;

    // path of parameters
    p_est[k]  = p;
    h2_est[k] = h2;
  }

  return List::create(
    _["beta_est"]       = avg_beta     / num_iter,
    _["postp_est"]      = avg_postp    / num_iter,
    _["corr_est"]       = avg_beta_hat / num_iter,
    _["path_p_est"]     = p_est,
    _["path_h2_est"]    = h2_est);
}

/******************************************************************************/
