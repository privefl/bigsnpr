/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

using namespace Rcpp;

/******************************************************************************/

// the one for bed files is in bed-fun.cpp

// [[Rcpp::export]]
List prod_and_rowSumsSq2(Environment BM,
                         const IntegerVector& ind_row,
                         const IntegerVector& ind_col,
                         const NumericVector& center,
                         const NumericVector& scale,
                         const NumericMatrix& V) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, ind_row, ind_col, BM["code256"], 1);

  int n = macc.nrow();
  int m = macc.ncol();
  myassert_size(m, V.rows());
  myassert_size(m, center.size());
  myassert_size(m, scale.size());
  int K = V.cols();

  NumericMatrix XV(n, K);
  NumericVector rowSumsSq(n);

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      double x = (macc(i, j) - center[j]) / scale[j];
      rowSumsSq[i] += x*x;
      for (int k = 0; k < K; k++) {
        XV(i, k) += x * V(j, k);
      }
    }
  }

  return List::create(XV, rowSumsSq);
}

/******************************************************************************/
