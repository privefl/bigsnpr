// [[Rcpp::depends(bigmemory, BH, RcppArmadillo)]]
#include <RcppArmadillo.h> // Sys.setenv("PKG_LIBS" = "-llapack")
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;


arma::mat getXtW(const arma::mat& covar, const arma::vec& w) {
  return trans(covar.each_col() % w);
}

// [[Rcpp::export]]
ListOf<SEXP> wcrossprod(SEXP pBigMat,
                        arma::mat& covar,
                        const arma::vec& y,
                        const arma::vec& z0,
                        const arma::vec& w0,
                        double tol,
                        int maxiter) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int n = xpMat->nrow();
  int m = xpMat->ncol();
  arma::mat tcovar, tmp;
  arma::vec p, w, z, betas_old, betas_new, Xb;
  double diff;
  int c;

  NumericVector res(m);
  NumericVector var(m);
  LogicalVector conv(m);

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      covar(i, 0) = macc[j][i];
    }
    z = z0;
    w = w0;
    tcovar = getXtW(covar, w);
    betas_new = solve(tcovar * covar, tcovar * z);
    c = 0;

    do {
      c++;
      betas_old = betas_new;

      Xb = covar * betas_old;
      p = 1 / (1 + exp(-Xb));
      w = p % (1 - p);
      z = Xb + (y - p) / w;

      tcovar = getXtW(covar, w);
      betas_new = solve(tcovar * covar, tcovar * z);

      diff = 2 * abs(betas_old(0) - betas_new(0))
        / (abs(betas_old(0)) + abs(betas_new(0)));
    } while (diff > tol && c < maxiter);

    res[j] = betas_new(0);
    tmp = inv(tcovar * covar);
    var[j] = tmp(0, 0);
    conv[j] = (c < maxiter);
  }

  return(List::create(_["betas"] = res,
                      _["std"] = sqrt(var),
                      _["conv"] = conv));
}
