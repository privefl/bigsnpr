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
List lassosum2(Environment corr,
               const NumericVector& beta_hat,
               double lambda,
               double s,
               double dfmax,
               int maxiter,
               double tol) {

  XPtr<SFBM> sfbm = corr["address"];

  int m = beta_hat.size();
  myassert_size(sfbm->nrow(), m);
  myassert_size(sfbm->ncol(), m);

  arma::vec curr_beta(m, arma::fill::zeros);

  double one_minus_s = 1 - s;
  if (one_minus_s == 0) dfmax = R_PosInf;

  int k = 0;
  for (; k < maxiter; k++) {

    double max_update = 0;
    double df = 0;

    for (int j = 0; j < m; j++) {

      double dotprod = sfbm->dot_col(j, curr_beta);
      double u_j = beta_hat[j] + one_minus_s * (curr_beta[j] - dotprod);
      double new_beta_j = soft_thres(u_j, lambda);
      if (new_beta_j != 0) df++;

      double shift = new_beta_j - curr_beta[j];
      double update = shift * shift;
      if (update > max_update) max_update = update;

      curr_beta[j] = new_beta_j;
    }

    if (max_update < tol || df > dfmax) {
      // Rcout << "Stopped after " << k + 1 << " iterations." << std::endl;
      break;
    }
  }

  return List::create(_["beta_est"] = curr_beta, _["num_iter"] = k + 1);
}

/******************************************************************************/
