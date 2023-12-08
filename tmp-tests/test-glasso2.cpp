// Based on algo described in DOI: 10.1371/journal.pone.0014147

#include <RcppArmadillo.h>
using namespace Rcpp;

inline double soft_thres(double z, double l1, double l2) {
  if (z > 0) {
    double num = z - l1;
    return (num > 0) ? num / l2 : 0;
  } else {
    double num = z + l1;
    return (num < 0) ? num / l2 : 0;
  }
}

void inner_lasso(const arma::mat& mat,
                 arma::mat& W,
                 arma::mat& beta,
                 arma::vec& dotprods,
                 double lambda,
                 int m,
                 int i,
                 int maxiter,
                 double tol) {

  // arma::vec dotprods = W * beta.col(i);  // use parallelism to fasten
  dotprods.fill(0);
  for (int j = 0; j < m; j++) {
    if (beta(j, i) != 0) {  // use sparsity to fasten
      dotprods += W.col(j) * beta(j, i);
    }
  }

  double gap0 = std::inner_product(mat.begin_col(i), mat.end_col(i),
                                   mat.begin_col(i), 0.0);

  for (int k = 0; k < maxiter; k++) {

    bool conv_inner = true;
    double gap = 0;

    for (int j = 0; j < m; j++) {

      if (j != i) {

        double resid = mat(j, i) - dotprods[j];
        gap += resid * resid;
        double curr_beta = beta(j, i);
        double new_beta = soft_thres(resid + curr_beta * W(j, j), lambda, W(j, j));

        double shift = new_beta - curr_beta;
        if (shift != 0) {
          if (conv_inner && std::abs(shift) > tol) conv_inner = false;
          beta(j, i) = new_beta;
          dotprods += W.col(j) * shift;
        }
      }
    }

    if (gap > gap0) Rcpp::stop("Divergence!");
    if (conv_inner) break;
  }

  for (int j = 0; j < m; j++)
    if (j != i)
      W(i, j) = W(j, i) = dotprods[j];
}

// [[Rcpp::export]]
ListOf<NumericMatrix> glasso(const arma::mat& mat,
                             double lambda,
                             int maxiter_outer,
                             int maxiter_lasso,
                             double tol,
                             bool verbose) {

  int m = mat.n_cols;

  arma::mat W = mat + 0;
  W.diag() += lambda;

  arma::mat beta(m, m, arma::fill::zeros);
  arma::vec dotprods(m);

  for (int k = 0; k < maxiter_outer; k++) {

    Rcpp::checkUserInterrupt();
    if (verbose) Rcpp::Rcout << k + 1 << std::endl;

    arma::mat Wold = W + 0;

    for (int i = 0; i < m; i++) {
      inner_lasso(mat, W, beta, dotprods, lambda, m, i, maxiter_lasso, tol);
    }

    double max_diff = max(max(abs(W - Wold)));
    double diff2 = mean(mean(square(W - Wold)));
    if (verbose) Rcpp::Rcout << max_diff << " // " << diff2 << std::endl;

    if (max_diff < tol) break;
  }

  W.diag() -= lambda;
  return List::create(as<NumericMatrix>(wrap(W)),
                      as<NumericMatrix>(wrap(beta)));
}
