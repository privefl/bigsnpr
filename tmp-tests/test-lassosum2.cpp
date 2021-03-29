/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/utils.h>
#include <bigsparser/SFBM.h>

/******************************************************************************/

inline double soft_thres(double z, double l1) {
  if (z > 0) {
    double num = z - l1;
    return (num > 0) ? num : 0;
  } else {
    double num = z + l1;
    return (num < 0) ? num : 0;
  }
}

/******************************************************************************/

// [[Rcpp::export]]
arma::vec lassosum2(Environment corr,
                    const NumericVector& beta_hat,
                    const NumericVector& beta_init,
                    const IntegerVector& order,
                    double lambda,
                    double s,
                    int maxiter = 200,
                    double tol = 1e-8) {

  XPtr<SFBM> sfbm = corr["address"];

  int m = beta_hat.size();
  myassert_size(sfbm->nrow(), m);
  myassert_size(sfbm->ncol(), m);
  myassert_size(order.size(), m);
  myassert_size(beta_init.size(), m);

  arma::vec curr_beta(beta_init.begin(), m);

  double one_minus_s = 1 - s;

  for (int k = 0; k < maxiter; k++) {

    double max_update = 0;

    for (const int& j : order) {

      double dotprod = sfbm->dot_col(j, curr_beta);
      double u_j = beta_hat[j] + one_minus_s * (curr_beta[j] - dotprod);
      double new_beta_j = soft_thres(u_j, lambda);

      double shift = new_beta_j - curr_beta[j];
      double update = shift * shift;
      if (update > max_update) max_update = update;

      curr_beta[j] = new_beta_j;
    }

    if (max_update < tol) {
      Rcout << "Stopped after " << k + 1 << " iterations." << std::endl;
      break;
    }
  }

  return curr_beta;
}

/******************************************************************************/
