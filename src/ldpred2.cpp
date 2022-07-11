/******************************************************************************/

#include <bigsparser/SFBM.h>
#include <bigstatsr/utils.h>

/******************************************************************************/

// [[Rcpp::export]]
NumericVector ldpred2_gibbs_one(Environment corr,
                                const NumericVector& beta_hat,
                                const NumericVector& beta_init,
                                const IntegerVector& order,
                                const NumericVector& n_vec,
                                double h2,
                                double p,
                                bool sparse,
                                int burn_in,
                                int num_iter) {

  XPtr<SFBM> sfbm = corr["address"];

  int m = beta_hat.size();
  myassert_size(sfbm->nrow(), m);
  myassert_size(sfbm->ncol(), m);
  myassert_size(order.size(), m);
  myassert_size(beta_init.size(), m);
  myassert_size(n_vec.size(), m);

  NumericVector curr_beta = Rcpp::clone(beta_init);
  NumericVector dotprods  = sfbm->prod(curr_beta);
  NumericVector avg_beta(m);

  double h2_per_var = h2 / (m * p);
  double inv_odd_p = (1 - p) / p;
  double gap0 =
    std::inner_product(beta_hat.begin(), beta_hat.end(), beta_hat.begin(), 0.0);

  for (int k = -burn_in; k < num_iter; k++) {

    double gap = 0;

    for (const int& j : order) {

      // double dotprod = sfbm->dot_col(j, curr_beta);
      double resid = beta_hat[j] - dotprods[j];
      gap += resid * resid;
      double res_beta_hat_j = curr_beta[j] + resid;

      double C1 = h2_per_var * n_vec[j];
      double C2 = 1 / (1 + 1 / C1);
      double C3 = C2 * res_beta_hat_j;
      double C4 = C2 / n_vec[j];

      double post_p_j = 1 /
        (1 + inv_odd_p * ::sqrt(1 + C1) * ::exp(-C3 * C3 / C4 / 2));

      double diff = -curr_beta[j];
      if (sparse && (post_p_j < p)) {
        curr_beta[j] = 0;
      } else {
        curr_beta[j] = (post_p_j > ::unif_rand()) ? ::Rf_rnorm(C3, ::sqrt(C4)) : 0;
        if (k >= 0) avg_beta[j] += C3 * post_p_j;
      }
      diff += curr_beta[j];
      if (diff != 0) dotprods = sfbm->incr_mult_col(j, dotprods, diff);
    }

    if (gap > gap0) { avg_beta.fill(NA_REAL); return avg_beta; }
  }

  return avg_beta / num_iter;
}

/******************************************************************************/
