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

arma::vec ldpred2_gibbs_one(const arma::sp_mat& corr,
                            const NumericVector& betas_hat,
                            const NumericVector& betas_init,
                            const NumericVector& order,
                            const NumericVector& n_vec,
                            double h2,
                            double p,
                            bool sparse,
                            int burn_in,
                            int num_iter) {

  int m = betas_hat.size();
  arma::vec curr_betas(betas_init.begin(), m);
  arma::vec curr_post_means(m, arma::fill::zeros);
  arma::vec avg_betas(m, arma::fill::zeros);

  int num_iter_tot = burn_in + num_iter;
  for (int k = 0; k < num_iter_tot; k++) {

    for (const int& j : order) {

      curr_betas[j] = 0;
      double res_beta_hat_j = betas_hat[j] - arma::dot(corr.col(j), curr_betas);


      double C1 = h2 * n_vec[j] / (m * p);
      double C2 = 1 / (1 + 1 / C1);
      double C3 = C2 * res_beta_hat_j ;
      double C4 = ::sqrt(C2 / n_vec[j]);

      double post_p_j = 1 /
        (1 + (1 - p) / p * ::sqrt(1 + C1) * ::exp(-square(C3 / C4) / 2));

      if (sparse && (post_p_j < p)) {
        curr_betas[j] = curr_post_means[j] = 0;
      } else {
        curr_post_means[j] = C2 * post_p_j * res_beta_hat_j;
        curr_betas[j] = (post_p_j > ::unif_rand()) ? (C3 + ::norm_rand() * C4) : 0;
      }
    }

    if (k >= burn_in) avg_betas += curr_post_means;
  }

  return avg_betas / num_iter;
}

/******************************************************************************/

// [[Rcpp::export]]
arma::mat ldpred2_gibbs(const arma::sp_mat& corr,
                        const NumericVector& betas_hat,
                        const NumericVector& betas_init,
                        const NumericVector& order,
                        const NumericVector& n_vec,
                        const NumericVector& h2,
                        const NumericVector& p,
                        const LogicalVector& sparse,
                        int burn_in,
                        int num_iter,
                        int ncores) {

  int m = betas_hat.size();
  myassert_size(corr.n_rows, m);
  myassert_size(corr.n_cols, m);
  myassert_size(betas_init.size(), m);
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
      corr, betas_hat, betas_init, order, n_vec,
      h2[k], p[k], sparse[k], burn_in, num_iter);

    #pragma omp critical
    res.col(k) = res_k;
  }

  return res;
}

/******************************************************************************/
