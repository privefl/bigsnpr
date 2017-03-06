// /******************************************************************************/
//
// #include "bigstatsr.h"
//
// /******************************************************************************/
//
// template <typename T>
// NumericMatrix linRegPcadapt(MatrixAccessor<T> macc,
//                             arma::mat &U,
//                             const IntegerVector &rowInd,
//                             const NumericVector &center,
//                             const NumericVector &scale) {
//   int n = rowInd.size();
//   int m = macc.ncol();
//   int K = U.n_cols;
//
//   arma::mat res(K, m);
//   arma::vec y(n), betas(K), eps(n);
//   int i, j;
//
//   // indices begin at 1 in R and 0 in C++
//   IntegerVector trains = rowInd - 1;
//
//   for (j = 0; j < m; j++) {
//     for (i = 0; i < n; i++) { // j-th SNP (centered)
//       y[i] = (macc[j][trains[i]] - center[j]) / scale[j];
//     }
//
//     betas = U.t() * y;
//     eps = y - U * betas;
//     res.col(j) = betas / sqrt(dot(eps, eps) / (n - K));
//   }
//
//   return as<NumericMatrix>(wrap(res.t()));
// }
//
// /******************************************************************************/
//
// // Dispatch function for linRegPcadapt
// //'
// //' T-scores used in pcadapt
// //'
// //' Compute matrix of t-scores (SNPs x scores) used in pcadapt.
// //'
// //' @param xpMat Slot `address` of a `big.matrix` object.
// //' @param U Matrix of left singular vectors (from partial SVD).
// //' @param rowInd Vector of row indices of the `big.matrix` that are used.
// //'
// //' @return A matrix of t-scores where rows correspond to each SNP and
// //' columns correspond to each left singular vector.
// //'
// //' @references Keurcien Luu and Michael Blum (2017).
// //' pcadapt: Fast Principal Component Analysis for Outlier Detection.
// //' R package version 3.0.4. https://CRAN.R-project.org/package=pcadapt.
// //'
// //' @export
// //' @keywords internal
// // [[Rcpp::export]]
// NumericMatrix linRegPcadapt(XPtr<BigMatrix> xpMat,
//                             arma::mat &U,
//                             const IntegerVector &rowInd,
//                             const NumericVector &center,
//                             const NumericVector &scale) {
//   myassert(U.n_rows == rowInd.size(), ERROR_DIM);
//
//   switch(xpMat->matrix_type()) {
//   case 1:
//     return linRegPcadapt(MatrixAccessor<char>(*xpMat),   U, rowInd,
//                          center, scale);
//   case 2:
//     return linRegPcadapt(MatrixAccessor<short>(*xpMat),  U, rowInd,
//                          center, scale);
//   case 4:
//     return linRegPcadapt(MatrixAccessor<int>(*xpMat),    U, rowInd,
//                          center, scale);
//   case 6:
//     return linRegPcadapt(MatrixAccessor<float>(*xpMat),  U, rowInd,
//                          center, scale);
//   case 8:
//     return linRegPcadapt(MatrixAccessor<double>(*xpMat), U, rowInd,
//                          center, scale);
//   default:
//     throw Rcpp::exception(ERROR_TYPE);
//   }
// }
//
// /******************************************************************************/
