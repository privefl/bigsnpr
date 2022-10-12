/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigsparser/SFBM.h>
#include <bigstatsr/utils.h>
#include "optim-MLE-alpha.h"

/******************************************************************************/

const double MIN_P  = 1e-5;
const double MIN_H2 = 1e-3;

/******************************************************************************/

double dotprod2(const NumericVector& X, const NumericVector& Y) {
  return std::inner_product(X.begin(), X.end(), Y.begin(), 0.0);
}

/******************************************************************************/

// [[Rcpp::export]]
arma::vec& MLE_alpha(arma::vec& par,
                     const std::vector<int>& ind_causal,
                     const NumericVector& log_var,
                     const NumericVector& curr_beta,
                     const NumericVector& alpha_bounds,
                     bool boot = false,
                     bool verbose = false) {

  MLE mle(ind_causal, log_var, curr_beta, boot);
  Roptim<MLE> opt("L-BFGS-B");
  opt.set_lower({alpha_bounds[0], par[1] / 2});
  opt.set_upper({alpha_bounds[1], par[1] * 2});
  opt.set_hessian(false);
  opt.control.trace = verbose;

  if (verbose) {
    arma::vec grad1(2), grad2(2);
    mle.Gradient(par, grad1);
    mle.ApproximateGradient(par, grad2);

    Rcout << "-------------------------" << std::endl;
    Rcout << "Gradient checking" << std::endl;
    grad1.t().print("analytic:");
    grad2.t().print("approximate:");
    Rcout << "-------------------------" << std::endl;
  }

  opt.minimize(mle, par);

  if (verbose) {
    Rcout << "-------------------------" << std::endl;
    opt.print();
    Rcout << "-------------------------" << std::endl;
  }

  return par;
}

/******************************************************************************/

// [[Rcpp::export]]
List ldpred2_gibbs_auto(Environment corr,
                        const NumericVector& beta_hat,
                        const NumericVector& beta_init,
                        const IntegerVector& order,
                        const NumericVector& n_vec,
                        const NumericVector& log_var,
                        double p_init,
                        double h2_init,
                        int burn_in,
                        int num_iter,
                        int report_step,
                        bool no_jump_sign,
                        double shrink_corr,
                        const NumericVector& alpha_bounds,
                        double mean_ld = 1,
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

  arma::sp_mat sample_beta(m, num_iter / report_step);
  int ind_report = 0, next_k_reported = burn_in + report_step - 1;

  int num_iter_tot = burn_in + num_iter;
  NumericVector p_est(num_iter_tot), h2_est(num_iter_tot), alpha_est(num_iter_tot);

  double cur_h2_est = shrink_corr * dotprod2(curr_beta, dotprods) +
    (1 - shrink_corr) * dotprod2(curr_beta, curr_beta);
  double p = std::max(p_init, MIN_P), h2 = std::max(h2_init, MIN_H2);
  arma::vec par_mle = {0, h2 / (m * p)};  // (alpha + 1) and sigma2 [init]
  std::vector<int> ind_causal;

  for (int k = 0; k < num_iter_tot; k++) {

    double inv_odd_p = (1 - p) / p;
    double alpha_plus_one = par_mle[0], sigma2 = par_mle[1];

    ind_causal.clear();

    for (const int& j : order) {

      double dotprod = dotprods[j];  // sfbm->dot_col(j, curr_beta);
      double res_beta_hat_j = beta_hat[j] + shrink_corr * (curr_beta[j] - dotprod);

      double scale_freq = ::exp(alpha_plus_one * log_var[j]);
      double C1 = scale_freq * sigma2 * n_vec[j];
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
          // if (k >= burn_in) {  // not sure if I should undo this
          //   avg_postp[j] -= postp;
          //   avg_beta[j]  -= C3 * postp;
          // }
        } else {
          curr_beta[j] = samp_beta;
          diff += samp_beta;
          ind_causal.push_back(j);
        }

      } else {
        curr_beta[j] = 0;
      }

      if (diff != 0) {
        cur_h2_est += diff * (2 * dotprod_shrunk + diff);
        dotprods = sfbm->incr_mult_col(j, dotprods, diff);
      }
    }

    int nb_causal = ind_causal.size();
    p = std::max(::Rf_rbeta(1 + nb_causal / mean_ld, 1 + (m - nb_causal) / mean_ld), MIN_P);
    h2 = std::max(cur_h2_est, MIN_H2);
    par_mle = MLE_alpha(par_mle, ind_causal, log_var, curr_beta, alpha_bounds, true);
    if (verbose) Rcout <<
      k + 1 << ": " << p << " // " << h2 << " // " <<  par_mle[0] - 1 << std::endl;

    // path of parameters
    p_est[k]  = p;
    h2_est[k] = h2;
    alpha_est[k] = par_mle[0] - 1;

    // store some sampling betas
    if (k == next_k_reported) {
      for (const int& i : ind_causal) {
        sample_beta(i, ind_report) = curr_beta[i];
      }
      ind_report++;
      next_k_reported += report_step;
    }
  }

  return List::create(
    _["beta_est"]       = avg_beta     / num_iter,
    _["postp_est"]      = avg_postp    / num_iter,
    _["corr_est"]       = avg_beta_hat / num_iter,
    _["sample_beta"]    = sample_beta,
    _["path_p_est"]     = p_est,
    _["path_h2_est"]    = h2_est,
    _["path_alpha_est"] = alpha_est);
}

/******************************************************************************/
