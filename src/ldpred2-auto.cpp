/******************************************************************************/

// [[Rcpp::plugins(cpp11)]]
#define ARMA_64BIT_WORD

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/utils.h>

using namespace Rcpp;

/******************************************************************************/

double h2_boot(const arma::vec& curr_betas) {

  double h2_est = 0;
  int m = curr_betas.size();
  double m2 = m - 1e-8;

  for (int k = 0; k < m; k++) {
    int i = ::unif_rand() * m2;
    double beta = curr_betas[i];
    h2_est += beta * beta;
  }

  return h2_est;
}

double p_boot(const std::vector<double>& post_p) {

  double p_est = 0;
  int m = post_p.size();
  double m2 = m - 1e-8;

  for (int k = 0; k < m; k++) {
    int i = ::unif_rand() * m2;
    double p = post_p[i];
    p_est += p;
  }

  return p_est / m;
}

/******************************************************************************/

List ldpred2_gibbs_auto_one(const arma::sp_mat& corr,
                            const NumericVector& betas_hat,
                            const NumericVector& betas_init,
                            const IntegerVector& order,
                            const NumericVector& n_vec,
                            double h2_max,
                            double p_init,
                            int burn_in,
                            int num_iter) {

  int m = betas_hat.size();
  std::vector<double> post_p(m, p_init);
  arma::vec curr_betas(m, arma::fill::zeros);
  std::copy(betas_init.begin(), betas_init.end(), curr_betas.begin());
  arma::vec curr_post_means(m, arma::fill::zeros);
  arma::vec avg_betas(m, arma::fill::zeros);

  int num_iter_tot = burn_in + num_iter;
  std::vector<double> p_est(num_iter_tot), h2_est(num_iter_tot);

  double p = p_init, h2 = h2_max, alpha = 1;
  double avg_p = 0, avg_h2 = 0;

  for (int k = 0; k < num_iter_tot; k++) {

    for (const int& j : order) {

      curr_betas[j] = 0;
      double res_beta_hat_j = betas_hat[j] - arma::dot(corr.col(j), curr_betas);

      double L = h2 * n_vec[j] / (m * p);
      double C2 = 1 / (1 + 1 / L);

      post_p[j] = 1 / (1 + (1 - p) / p * ::sqrt(1 + L) *
        ::exp(-C2 * n_vec[j] / 2 * res_beta_hat_j * res_beta_hat_j));

      curr_post_means[j] = C2 * post_p[j] * res_beta_hat_j;
      curr_betas[j] = ((alpha * post_p[j]) > ::unif_rand())
        ? C2 * res_beta_hat_j + ::norm_rand() * ::sqrt(C2 / n_vec[j]) : 0;
    }

    p     = std::max(1e-5, p_boot(post_p));
    h2    = std::max(1e-4, h2_boot(curr_betas));
    alpha = std::min(1.0,  h2_max / h2);

    if (k >= burn_in) {
      avg_betas += curr_post_means;
      avg_p += p;
      avg_h2 += h2;
    }
    Rcout << k + 1 << ": " << p << " // " << h2 << std::endl;
    p_est[k] = p;
    h2_est[k] = h2;
  }

  Rcout << (avg_p / num_iter) << " // " << (avg_h2 / num_iter) << std::endl;

  return List::create(
    _["beta_est"]   = avg_betas / num_iter,
    _["p_est"]      = avg_p     / num_iter,
    _["h2_est"]     = avg_h2    / num_iter,
    _["vec_p_est"]  = wrap(p_est),
    _["vec_h2_est"] = wrap(h2_est));
}

/******************************************************************************/

// [[Rcpp::export]]
List ldpred2_gibbs_auto(const arma::sp_mat& corr,
                        const NumericVector& betas_hat,
                        const NumericVector& betas_init,
                        const IntegerVector& order,
                        const NumericVector& n_vec,
                        const NumericVector& h2_max,
                        const NumericVector& p_init,
                        int burn_in,
                        int num_iter,
                        int ncores = 1) {

  myassert_size(h2_max.size(), p_init.size());
  myassert_size(corr.n_cols, betas_hat.size());
  myassert_size(corr.n_cols, order.size());
  myassert_size(corr.n_cols, n_vec.size());

  int K = h2_max.size();
  List res(K);

  #pragma omp parallel for schedule(dynamic, 1) num_threads(ncores)
  for (int k = 0; k < K; k++) {
    List res_k = ldpred2_gibbs_auto_one(
      corr, betas_hat, betas_init, order, n_vec, h2_max[k], p_init[k], burn_in, num_iter);
    #pragma omp critical
    res[k] = res_k;
  }

  return res;
}

/******************************************************************************/
