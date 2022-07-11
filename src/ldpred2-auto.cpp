/******************************************************************************/

#include <bigsparser/SFBM.h>
#include <bigstatsr/utils.h>

/******************************************************************************/

const double MIN_P  = 1e-5;
const double MIN_H2 = 1e-4;

/******************************************************************************/

double dotprod2(const NumericVector& X, const NumericVector& Y) {
  return std::inner_product(X.begin(), X.end(), Y.begin(), 0.0);
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
                        int report_step,
                        bool allow_jump_sign,
                        double shrink_corr,
                        bool verbose = false) {

  XPtr<SFBM> sfbm = corr["address"];

  int m = beta_hat.size();
  myassert_size(sfbm->nrow(), m);
  myassert_size(sfbm->ncol(), m);
  myassert_size(order.size(), m);
  myassert_size(beta_init.size(), m);
  myassert_size(n_vec.size(), m);

  NumericVector curr_beta = Rcpp::clone(beta_init);
  NumericVector dotprods  = sfbm->prod(curr_beta);
  NumericVector avg_beta(m), avg_postp(m), avg_beta_hat(m);

  NumericMatrix sample_beta(m, num_iter / report_step);
  int ind_report = 0, next_k_reported = burn_in - 1 + report_step;

  int num_iter_tot = burn_in + num_iter;
  NumericVector p_est(num_iter_tot), h2_est(num_iter_tot);

  double cur_h2_est = shrink_corr * dotprod2(curr_beta, dotprods) +
    (1 - shrink_corr) * dotprod2(curr_beta, curr_beta);
  double p = p_init, h2 = h2_init, avg_p = 0, avg_h2 = 0;
  bool no_jump_sign = !allow_jump_sign;

  for (int k = 0; k < num_iter_tot; k++) {

    int nb_causal = 0;
    double h2_per_var = h2 / (m * p);
    double inv_odd_p = (1 - p) / p;

    for (const int& j : order) {

      double dotprod = dotprods[j];  // sfbm->dot_col(j, curr_beta);
      double res_beta_hat_j = beta_hat[j] + shrink_corr * (curr_beta[j] - dotprod);

      double C1 = h2_per_var * n_vec[j];
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
          if (k >= burn_in) {  // not sure if I should undo this
            avg_postp[j] -= postp;
            avg_beta[j]  -= C3 * postp;
          }
        } else {
          curr_beta[j] = samp_beta;
          diff += samp_beta;
          nb_causal++;
        }

      } else {
        curr_beta[j] = 0;
      }

      if (diff != 0) {
        cur_h2_est += diff * (2 * dotprod_shrunk + diff);
        dotprods = sfbm->incr_mult_col(j, dotprods, diff);
      }
    }

    p = std::max(::Rf_rbeta(1 + nb_causal, 1 + m - nb_causal), MIN_P);
    h2 = std::max(cur_h2_est, MIN_H2);
    if (verbose) Rcout << k + 1 << ": " << p << " // " << h2 << std::endl;

    if (k >= burn_in) {
      avg_p  += p;
      avg_h2 += h2;
      if (k == next_k_reported) {
        sample_beta(_, ind_report++) = curr_beta;
        next_k_reported += report_step;
      }
    }
    p_est[k]  = p;
    h2_est[k] = h2;
  }

  double est_p  = avg_p  / num_iter;
  double est_h2 = avg_h2 / num_iter;
  if (verbose) Rcout << "Overall: " << est_p << " // " << est_h2 << std::endl;

  return List::create(
    _["beta_est"]    = avg_beta     / num_iter,
    _["postp_est"]   = avg_postp    / num_iter,
    _["corr_est"]    = avg_beta_hat / num_iter,
    _["sample_beta"] = sample_beta,
    _["p_est"]       = est_p,
    _["h2_est"]      = est_h2,
    _["path_p_est"]  = p_est,
    _["path_h2_est"] = h2_est,
    _["h2_init"]     = h2_init,
    _["p_init"]      = p_init);
}

/******************************************************************************/
