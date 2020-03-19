/******************************************************************************/

// [[Rcpp::plugins(cpp11)]]
#define ARMA_64BIT_WORD

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/utils.h>

using namespace Rcpp;

/******************************************************************************/

arma::vec ldpred2_gibbs_one(const arma::sp_mat& corr,
                            const NumericVector& betas_hat,
                            const IntegerVector& order,
                            const NumericVector& n_vec,
                            double h2,
                            double p,
                            int burn_in,
                            int num_iter,
                            double sparse) {

  int m = betas_hat.size();
  arma::vec curr_betas(m, arma::fill::zeros);
  arma::vec curr_post_means(m, arma::fill::zeros);
  arma::vec avg_betas(m, arma::fill::zeros);

  int num_iter_tot = burn_in + num_iter;
  for (int k = 0; k < num_iter_tot; k++) {

    for (const int& j : order) {

      curr_betas[j] = 0;
      double res_beta_hat_j = betas_hat[j] - arma::dot(corr.col(j), curr_betas);

      double L = h2 * n_vec[j] / (m * p);
      double C2 = 1 / (1 + 1 / L);

      double post_p_j = 1 / (1 + (1 - p) / p * ::sqrt(1 + L) *
        ::exp(-C2 * n_vec[j] / 2 * res_beta_hat_j * res_beta_hat_j));

      if (sparse && (post_p_j < p)) {
        curr_betas[j] = curr_post_means[j] = 0;
      } else {
        curr_post_means[j] = C2 * post_p_j * res_beta_hat_j;
        curr_betas[j] = (post_p_j > ::unif_rand())
          ? C2 * res_beta_hat_j + ::norm_rand() * ::sqrt(C2 / n_vec[j]) : 0;
      }
    }

    if (k >= burn_in) avg_betas += curr_post_means;
  }

  return avg_betas / num_iter;
}

/******************************************************************************/

// [[Rcpp::export]]
List ldpred2_gibbs(const arma::sp_mat& corr,
                   const NumericVector& betas_hat,
                   const IntegerVector& order,
                   const NumericVector& n_vec,
                   const NumericVector& h2,
                   const NumericVector& p,
                   int burn_in,
                   int num_iter,
                   double sparse,
                   int ncores) {

  myassert_size(h2.size(), p.size());
  myassert_size(corr.n_cols, betas_hat.size());
  myassert_size(corr.n_cols, order.size());
  myassert_size(corr.n_cols, n_vec.size());

  int K = h2.size();
  List res(K);

  #pragma omp parallel for schedule(dynamic, 1) num_threads(ncores)
  for (int k = 0; k < K; k++) {
    arma::vec res_k = ldpred2_gibbs_one(
      corr, betas_hat, order, n_vec, h2[k], p[k], burn_in, num_iter, sparse);
    #pragma omp critical
    res[k] = wrap(res_k);
  }

  return res;
}

/******************************************************************************/
