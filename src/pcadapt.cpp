/******************************************************************************/

#include "bigsnpr.h"

/******************************************************************************/

//' T-scores used in pcadapt
//'
//' Compute matrix of t-scores (SNPs x scores) used in pcadapt.
//'
//' @param xpMat Slot `address` of a `big.matrix` object.
//' @param U Matrix of left singular vectors (from partial SVD).
//' @param rowInd Vector of row indices of the `big.matrix` that are used.
//'
//' @return A matrix of t-scores where rows correspond to each SNP and
//' columns correspond to each left singular vector.
//'
//' @references Keurcien Luu and Michael Blum (2017).
//' pcadapt: Fast Principal Component Analysis for Outlier Detection.
//' R package version 3.0.4. https://CRAN.R-project.org/package=pcadapt.
//'
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix linRegPcadapt(const S4& BM,
                            arma::mat& U,
                            const IntegerVector& rowInd) {



  XPtr<BigMatrix> xpMat = BM.slot("address");
  int m = xpMat->ncol();
  RawSubMatAcc macc(*xpMat, rowInd-1, seq_len(m)-1, BM.slot("code"));
  int n = macc.nrow();
  int K = U.n_cols;
  myassert((int)U.n_rows == n, ERROR_DIM);

  arma::mat res(K, m);
  arma::vec y(n), betas(K), eps(n);
  double sum_y;
  int i, j;

  for (j = 0; j < m; j++) {
    sum_y = 0;
    for (i = 0; i < n; i++) { // j-th SNP (centered)
      y[i] = macc(i, j);
      sum_y += y[i];
    }
    y -= sum_y / n;

    betas = U.t() * y;
    eps = y - U * betas;
    res.col(j) = betas / sqrt(dot(eps, eps) / (n - K));
  }

  return as<NumericMatrix>(wrap(res.t()));
}

/******************************************************************************/
