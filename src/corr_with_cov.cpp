// /******************************************************************************/
//
// #include <bigstatsr/BMCodeAcc.h>
// #include "bed-acc.h"
//
// using namespace Rcpp;
//
// /******************************************************************************/
//
// template <class C>
// List corMat0(C macc,
//              const NumericMatrix& Z,
//              double size,
//              const NumericVector& thr,
//              const NumericVector& pos,
//              int ncores) {
//
//   int n = macc.nrow();
//   int m = macc.ncol();
//   int K = Z.nrow();
//   int df = n - K;
//
//   List res(m);
//   std::vector<double> sqrt_var(m, NA_REAL);
//
//   int chunk_size = ceil(m / (10.0 * ncores));
//
//   #pragma omp parallel num_threads(ncores)
//   {
//     std::vector<int>    ind; ind.reserve(m);
//     std::vector<double> val; val.reserve(m);
//
//     #pragma omp for schedule(dynamic, chunk_size)
//     for (int j0 = 0; j0 < m; j0++) {
//
//       ind.clear();
//       val.clear();
//
//       // pre-computation
//       double xxSum = 0;
//       for (int i = 0; i < n; i++) xxSum += macc(i, j0) * macc(i, j0);
//       for (int k = 0; k < K; k++) xxSum -=    Z(k, j0) *    Z(k, j0);
//
//       // main computation
//       double pos_min = pos[j0] - size;
//       for (int j = j0 - 1; (j >= 0) && (pos[j] >= pos_min); j--) {
//
//         double xySum = 0, yySum = 0;
//
//         for (int i = 0; i < n; i++) {
//           double y = macc(i, j);
//           xySum += macc(i, j0) * y;
//           yySum += y * y;
//         }
//         for (int k = 0; k < K; k++) {
//           double y = Z(k, j);
//           xySum -= Z(k, j0) * y;
//           yySum -= y * y:
//         }
//
//         double r = xySum / ::sqrt(xxSum * yySum);
//
//         if (ISNAN(r) || std::abs(r) > thr[df]) {
//           ind.push_back(j + 1);
//           val.push_back(r);
//         }
//       }
//
//       #pragma omp critical
//       res[j0] = List::create(_["i"] = wrap(ind), _["x"] = wrap(val));
//     }
//   }
//
//   return res;
// }
//
// /******************************************************************************/
//
// // [[Rcpp::export]]
// List corMat(Environment obj,
//             const IntegerVector& rowInd,
//             const IntegerVector& colInd,
//             const NumericMatrix& Z,
//             double size,
//             const NumericVector& thr,
//             const NumericVector& pos,
//             int ncores) {
//
//   myassert_size(colInd.size(), pos.size());
//
//   if (obj.exists("code256")) {
//     XPtr<FBM> xpBM = obj["address"];
//     SubBMCode256Acc macc(xpBM, rowInd, colInd, obj["code256"], 1);
//     return corMat0(macc, size, thr, pos, ncores);
//   } else if (obj.exists("bedfile")) {
//     XPtr<bed> xp_bed = obj["address"];
//     bedAcc macc(xp_bed, rowInd, colInd, NA_REAL);
//     return corMat0(macc, size, thr, pos, ncores);
//   } else {
//     throw Rcpp::exception("Unknown object type.");
//   }
// }
//
// /******************************************************************************/
