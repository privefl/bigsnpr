// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/isna.hpp>

using namespace Rcpp;


/******************************************************************************/

// // [[Rcpp::export]]
// void rawToBigPart(const IntegerMatrix& source,
//                   SEXP pBigMat,
//                   int colOffset = 0) {
//   XPtr<BigMatrix> xpMat(pBigMat);
//   MatrixAccessor<char> macc(*xpMat);
//
//   int nrows = source.rows();
//   int ncols = source.cols();
//   int tmp;
//
//   for (int j = 0; j < ncols; j++) {
//     for (int i = 0; i < nrows; i++) {
//       tmp = source(i, j);
//       if (isna(tmp)) {
//         macc[j+colOffset][i] = NA_CHAR;
//       } else {
//         macc[j+colOffset][i] = tmp;
//       }
//     }
//   }
//
//   return;
// }

/******************************************************************************/
