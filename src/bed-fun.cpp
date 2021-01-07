/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include "bed-acc.h"

/******************************************************************************/

// [[Rcpp::export]]
List bed_colstats(Environment obj_bed,
                  const IntegerVector& ind_row,
                  const IntegerVector& ind_col,
                  int ncores) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAcc macc(xp_bed, ind_row, ind_col);
  int n = macc.nrow();
  int m = macc.ncol();

  NumericVector sumX(m), denoX(m);
  IntegerVector nb_nona_col(m);

  #pragma omp parallel for num_threads(ncores)
  for (int j = 0; j < m; j++) {
    double xSum = 0, xxSum = 0;
    int c = n;
    for (int i = 0; i < n; i++) {
      double x = macc(i, j);
      if (x != 3) {
        xSum  += x;
        xxSum += x * x;
      } else {
        c--;
      }
    }
    sumX[j]  = xSum;
    denoX[j] = xxSum - xSum * xSum / c;
    nb_nona_col[j] = c;
  }

  int n_bad = Rcpp::sum(2 * nb_nona_col < n);
  if (n_bad > 0) Rcpp::warning("%d variants have >50%% missing values.", n_bad);

  return List::create(_["sumX"]  = sumX,
                      _["denoX"] = denoX,
                      _["nb_nona_col"] = nb_nona_col);
}

/******************************************************************************/

// [[Rcpp::export]]
IntegerMatrix bed_col_counts_cpp(Environment obj_bed,
                                 const IntegerVector& ind_row,
                                 const IntegerVector& ind_col,
                                 int ncores) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAcc macc(xp_bed, ind_row, ind_col);
  size_t n = macc.nrow();
  size_t m = macc.ncol();

  IntegerMatrix res(4, m);

  #pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) {
    for (size_t i = 0; i < n; i++) (res(macc(i, j), j))++;
  }

  return res;
}

// [[Rcpp::export]]
arma::Mat<int> bed_row_counts_cpp(Environment obj_bed,
                                 const IntegerVector& ind_row,
                                 const IntegerVector& ind_col,
                                 int ncores) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAcc macc(xp_bed, ind_row, ind_col);
  size_t n = macc.nrow();
  size_t m = macc.ncol();

  arma::Mat<int> res(4, n, arma::fill::zeros);

  #pragma omp parallel num_threads(ncores)
  {
    arma::Mat<int> res_local(4, n, arma::fill::zeros);

    #pragma omp for
    for (size_t j = 0; j < m; j++) {
      for (size_t i = 0; i < n; i++) (res_local(macc(i, j), i))++;
    }

    #pragma omp critical
    res += res_local;
  }

  return res;
}

/******************************************************************************/

// [[Rcpp::export]]
NumericMatrix read_bed_scaled(Environment obj_bed,
                              const IntegerVector& ind_row,
                              const IntegerVector& ind_col,
                              const NumericVector& center,
                              const NumericVector& scale) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAccScaled macc_bed(xp_bed, ind_row, ind_col, center, scale);

  size_t n = macc_bed.nrow();
  size_t m = macc_bed.ncol();

  NumericMatrix res(n, m);

  for (size_t j = 0; j < m; j++)
    for (size_t i = 0; i < n; i++)
      res(i, j) = macc_bed(i, j);

  return res;
}

/******************************************************************************/

// [[Rcpp::export]]
List prod_and_rowSumsSq(Environment obj_bed,
                        const IntegerVector& ind_row,
                        const IntegerVector& ind_col,
                        const NumericVector& center,
                        const NumericVector& scale,
                        const NumericMatrix& V) {

  XPtr<bed> xp_bed = obj_bed["address"];
  bedAccScaled macc_bed(xp_bed, ind_row, ind_col, center, scale);

  size_t n = macc_bed.nrow();
  size_t m = macc_bed.ncol();
  myassert_size(m, V.rows());
  size_t K = V.cols();
  size_t i, j, k;

  NumericMatrix XV(n, K);
  NumericVector rowSumsSq(n);

  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      double x = macc_bed(i, j);
      rowSumsSq[i] += x*x;
      for (k = 0; k < K; k++) {
        XV(i, k) += x * V(j, k);
      }
    }
  }

  return List::create(XV, rowSumsSq);
}

/******************************************************************************/

// // [[Rcpp::export]]
// NumericVector bed_wmean(Environment obj_bed,
//                         const IntegerVector& ind_row,
//                         const IntegerVector& ind_col,
//                         const NumericVector& w) {
//
//   myassert_size(w.size(), ind_row.size());
//
//   XPtr<bed> xp_bed = obj_bed["address"];
//   bedAcc macc(xp_bed, ind_row, ind_col);
//   size_t n = macc.nrow();
//   size_t m = macc.ncol();
//
//   NumericVector wmean(m);
//   double x, xwSum, wSum, wSum0 = Rcpp::sum(w);
//   int c;
//
//   for (size_t j = 0; j < m; j++) {
//     xwSum = 0;
//     wSum = wSum0;
//     c = n;
//     for (size_t i = 0; i < n; i++) {
//       x = macc(i, j);
//       if (x != 3) {
//         xwSum += x * w[i];
//       } else {
//         c--;
//         wSum -= w[i];
//       }
//     }
//     wmean[j] = xwSum / wSum;
//   }
//
//   return wmean;
// }

/******************************************************************************/

// // [[Rcpp::export]]
// List bed_corNA(Environment obj_bed,
//                const IntegerVector& ind_row,
//                const IntegerVector& ind_col,
//                const NumericMatrix& U) {
//
//   XPtr<bed> xp_bed = obj_bed["address"];
//   bedAcc macc(xp_bed, ind_row, ind_col);
//   size_t n = macc.nrow();
//   size_t m = macc.ncol();
//   myassert_size(U.nrow(), ind_row.size());
//
//   size_t K = U.ncol();
//
//   NumericMatrix res(K, m);
//   IntegerVector nb_na(m);
//
//   for (size_t j = 0; j < m; j++) {
//     for (size_t i = 0; i < n; i++) {
//       if (macc(i, j) == 3) {
//         nb_na[j]++;
//         for (size_t k = 0; k < K; k++) res(k, j) += U(i, k);
//       }
//     }
//   }
//
//   return List::create(_["prod"]  = res,
//                       _["nb_na"] = nb_na);
// }

/******************************************************************************/
