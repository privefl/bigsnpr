/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/utils.h>
#include <bigsparser/SFBM.h>

/******************************************************************************/

inline double soft_thres(double z, double l1, double one_plus_l2) {
  if (z > 0) {
    double num = z - l1;
    return (num > 0) ? num / one_plus_l2 : 0;
  } else {
    double num = z + l1;
    return (num < 0) ? num / one_plus_l2 : 0;
  }
}

/******************************************************************************/

// [[Rcpp::export]]
List lassosum2(Environment corr,
               const arma::vec& beta_hat,
               double lambda,
               double delta,
               double dfmax,
               int maxiter,
               double tol) {

  XPtr<SFBM> sfbm = corr["address"];

  int m = beta_hat.size();
  myassert_size(sfbm->nrow(), m);
  myassert_size(sfbm->ncol(), m);

  arma::vec curr_beta(m, arma::fill::zeros), dotprods(m, arma::fill::zeros);

  double one_plus_delta = 1 + delta;
  double gap0 = arma::dot(beta_hat, beta_hat);

  int k = 0;
  for (; k < maxiter; k++) {

    bool conv = true;
    double df = 0;
    double gap = 0;

    for (int j = 0; j < m; j++) {

      double resid = beta_hat[j] - dotprods[j];
      gap += resid * resid;
      double u_j = curr_beta[j] + resid;
      double new_beta_j = soft_thres(u_j, lambda, one_plus_delta);
      if (new_beta_j != 0) df++;

      double shift = new_beta_j - curr_beta[j];
      if (shift != 0) {
        if (conv && std::abs(shift) > tol) conv = false;
        curr_beta[j] = new_beta_j;
        dotprods = sfbm->incr_mult_col(j, dotprods, shift);
      }
    }

    if (gap > gap0) { curr_beta.fill(NA_REAL); break; }
    if (conv || df > dfmax) break;
  }

  return List::create(_["beta_est"] = curr_beta, _["num_iter"] = k + 1);
}

/******************************************************************************/
