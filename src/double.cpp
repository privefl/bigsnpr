/******************************************************************************/

#include "bigsnpr.h"

/******************************************************************************/

// [[Rcpp::export]]
void doubleBM(const S4& BM, XPtr<BigMatrix> xpMat2) {

  XPtr<BigMatrix> xpMat = BM.slot("address");
  int n = xpMat->nrow();
  int m = xpMat->ncol();
  RawSubMatAcc macc(*xpMat, seq_len(n)-1, seq_len(m)-1, BM.slot("code"));

  MatrixAccessor<unsigned char> macc2(*xpMat2);

  int i, j, j2;
  double tmp;

  for (j = j2 = 0; j < m; j++, j2 += 2) {
    for (i = 0; i < n; i++) {
      tmp = macc(i, j);
      if (tmp == 0) {
        macc2[j2][i] = macc2[j2+1][i] = 0;
      } else if (tmp == 1) {
        macc2[j2][i] = 0;
        macc2[j2+1][i] = 1;
      } else if (tmp == 2) {
        macc2[j2][i] = macc2[j2+1][i] = 1;
      } else {
        throw Rcpp::exception("Your big.matrix should have only Os, 1s or 2s");
      }
    }
  }
}

/******************************************************************************/

// [[Rcpp::export]]
void tripleBM(const S4& BM, XPtr<BigMatrix> xpMat2) {

  XPtr<BigMatrix> xpMat = BM.slot("address");
  int n = xpMat->nrow();
  int m = xpMat->ncol();
  RawSubMatAcc macc(*xpMat, seq_len(n)-1, seq_len(m)-1, BM.slot("code"));

  MatrixAccessor<unsigned char> macc2(*xpMat2);

  int i, j, j2;
  double tmp;

  for (j = j2 = 0; j < m; j++, j2 += 3) {
    for (i = 0; i < n; i++) {
      tmp = macc(i, j);
      if (tmp == 0) {
        macc2[j2][i] = macc2[j2+1][i] = macc2[j2+2][i] = 0;
      } else if (tmp == 1) {
        macc2[j2][i] = macc2[j2+2][i] = 1;
        macc2[j2+1][i] = 0;
      } else if (tmp == 2) {
        macc2[j2][i] = 2;
        macc2[j2+1][i] = macc2[j2+2][i] = 1;
      } else {
        throw Rcpp::exception("Your big.matrix should have only Os, 1s or 2s");
      }
    }
  }
}
