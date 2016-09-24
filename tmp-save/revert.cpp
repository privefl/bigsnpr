
// [[Rcpp::export]]
void revert(const SEXP pBigMat, SEXP pBigMat2, const LogicalVector& cond) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);
  XPtr<BigMatrix> xpMat2(pBigMat2);
  MatrixAccessor<char> macc2(*xpMat2);

  int n = xpMat->nrow();
  int m = xpMat->ncol();

  for (int j = 0; j < m; j++) {
    if (cond[j]) {
      for (int i = 0; i < n; i++) {
        macc2[j][i] = -macc[j][i] + 2;
      }
    } else {
      for (int i = 0; i < n; i++) {
        macc2[j][i] = macc[j][i];
      }
    }
  }

  return;
}
