/******************************************************************************/

#include <bigsparser/SFBM.h>
#include <bigstatsr/utils.h>

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
               const NumericVector& beta_hat,
               const NumericVector& lambda,
               const NumericVector& delta_plus_one,
               double dfmax,
               int maxiter,
               double tol) {

  XPtr<SFBM> sfbm = corr["address"];

  int m = beta_hat.size();
  myassert_size(sfbm->nrow(), m);
  myassert_size(sfbm->ncol(), m);

  NumericVector curr_beta(m), dotprods(m);

  double gap0 =
    std::inner_product(beta_hat.begin(), beta_hat.end(), beta_hat.begin(), 0.0);

  int k = 0;
  for (; k < maxiter; k++) {

    bool conv = true;
    double df = 0;
    double gap = 0;

    for (int j = 0; j < m; j++) {

      double resid = beta_hat[j] - dotprods[j];
      gap += resid * resid;
      double u_j = curr_beta[j] + resid;
      double new_beta_j = soft_thres(u_j, lambda[j], delta_plus_one[j]);
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
