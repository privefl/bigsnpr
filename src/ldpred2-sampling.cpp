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
arma::mat ldpred2_gibbs_one_sampling(Environment corr,
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

  arma::vec curr_beta(beta_init.begin(), m);
  arma::mat sample_beta(m, num_iter, arma::fill::zeros);

  int num_iter_tot = burn_in + num_iter;
  for (int k = 0; k < num_iter_tot; k++) {

    for (const int& j : order) {

      double dotprod = sfbm->dot_col(j, curr_beta);
      double res_beta_hat_j = beta_hat[j] + curr_beta[j] - dotprod;

      double C1 = h2 * n_vec[j] / (m * p);
      double C2 = 1 / (1 + 1 / C1);
      double C3 = C2 * res_beta_hat_j;
      double C4 = ::sqrt(C2 / n_vec[j]);

      double post_p_j = 1 /
        (1 + (1 - p) / p * ::sqrt(1 + C1) * ::exp(-square(C3 / C4) / 2));

      if (sparse && (post_p_j < p)) {
        curr_beta[j] = 0;
      } else {
        curr_beta[j] = (post_p_j > ::unif_rand()) ? (C3 + ::norm_rand() * C4) : 0;
        if (k >= burn_in) sample_beta(j, k - burn_in) = curr_beta[j];
      }
    }
  }

  return sample_beta;
}

/******************************************************************************/
