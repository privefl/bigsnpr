// [[Rcpp::depends(bigmemory, BH)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;


/******************************************************************************/

// // [[Rcpp::export]]
// NumericVector prs2(SEXP pBigMat,
//                    const NumericMatrix& odds,
//                    const IntegerVector& indTest,
//                    const IntegerVector& indCol) {
//   XPtr<BigMatrix> xpMat(pBigMat);
//   MatrixAccessor<char> macc(*xpMat);
//
//   int n = indTest.size();
//   int m = indCol.size();
//
//   // indices begin at 1 in R and 0 in C++
//   IntegerVector tests = indTest - 1;
//
//   NumericVector res(n);
//   int i, j, k;
//
//   for (j = 0; j < m; j++) {
//     k = indCol[j]-1;
//     for (i = 0; i < n; i++) {
//       res[i] += odds(macc[k][tests[i]], k);
//     }
//   }
//
//   return(res);
// }

/******************************************************************************/

// [[Rcpp::export]]
NumericVector prs1(XPtr<BigMatrix> xpMat,
                   const NumericVector& betas,
                   const IntegerVector& indTest,
                   const IntegerVector& indCol) {

  MatrixAccessor<char> macc(*xpMat);

  int n = indTest.size();
  int m = indCol.size();

  // indices begin at 1 in R and 0 in C++
  IntegerVector tests = indTest - 1;

  NumericVector res(n);
  int i, j, k;

  for (j = 0; j < m; j++) {
    k = indCol[j]-1;
    for (i = 0; i < n; i++) {
      res[i] += macc[k][tests[i]] * betas[k];
    }
  }

  return(res);
}

/******************************************************************************/
