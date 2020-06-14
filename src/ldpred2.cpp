/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/utils.h>
#include <bigsparser/SFBM.h>

/******************************************************************************/

inline double square(double x) {
  return x * x;
}

/******************************************************************************/

arma::vec ldpred2_gibbs_one(XPtr<SFBM> sfbm,
                            const NumericVector& beta_hat,
                            const NumericVector& beta_init,
                            const IntegerVector& order,
                            const NumericVector& n_vec,
                            double h2,
                            double p,
                            bool sparse,
                            int burn_in,
                            int num_iter) {

  int m = beta_hat.size();
  arma::vec curr_beta(beta_init.begin(), m);
  arma::vec avg_beta(m, arma::fill::zeros);

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
        if (k >= burn_in) avg_beta[j] += C3 * post_p_j;
      }
    }
  }

  return avg_beta / num_iter;
}

/******************************************************************************/

// [[Rcpp::export]]
arma::mat ldpred2_gibbs(Environment corr,
                        const NumericVector& beta_hat,
                        const NumericVector& beta_init,
                        const IntegerVector& order,
                        const NumericVector& n_vec,
                        const NumericVector& h2,
                        const NumericVector& p,
                        const LogicalVector& sparse,
                        int burn_in,
                        int num_iter,
                        int ncores) {

  XPtr<SFBM> sfbm = corr["address"];

  int m = beta_hat.size();
  myassert_size(sfbm->nrow(), m);
  myassert_size(sfbm->ncol(), m);
  myassert_size(order.size(), m);
  myassert_size(beta_init.size(), m);
  myassert_size(n_vec.size(), m);

  int K = p.size();
  myassert_size(h2.size(),     K);
  myassert_size(sparse.size(), K);

  arma::mat res(m, K);

  #pragma omp parallel for schedule(dynamic, 1) num_threads(ncores)
  for (int k = 0; k < K; k++) {

    // if (k % ncores == 0)
    //   Rcout << "Starting with params " << k + 1 << " / " << K << std::endl;

    arma::vec res_k = ldpred2_gibbs_one(
      sfbm, beta_hat, beta_init, order, n_vec,
      h2[k], p[k], sparse[k], burn_in, num_iter);

    #pragma omp critical
    res.col(k) = res_k;
  }

  return res;
}

/******************************************************************************/
