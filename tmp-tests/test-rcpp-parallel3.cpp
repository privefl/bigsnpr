// [[Rcpp::depends(rmio, bigstatsr)]]
#include <bigstatsr/BMCodeAcc.h>

// [[Rcpp::export]]
NumericVector parallelVectorSum(Environment BM) {

  XPtr<FBM> xpBM = BM["address"];
  std::size_t n = xpBM->nrow();
  std::size_t m = xpBM->ncol();
  SubBMCode256Acc macc(xpBM, seq_len(n), seq_len(m), BM["code256"], 1);

  NumericVector res(m);

  // std::size_t j;
  // #pragma omp parallel for private(j) schedule(static)
  // for (j = 0; j < m; j++) {
  //   double xySum = 0;
  //   for (std::size_t i = 0; i < n; i++) {
  //     xySum += macc(i, j) * macc(i, 0);
  //   }
  //   res[j] = xySum;
  // }

  for (std::size_t j = 0; j < m; j++) {
    double xySum = 0;
    std::size_t i;
    #pragma omp parallel for private(i) reduction(+:xySum) schedule(static, 1) num_threads(2)
    for (i = 0; i < n; i++) {
      xySum += macc(i, j) * macc(i, 0);
    }
    res[j] = xySum;
  }

  return res;
}


/*** R
library(bigsnpr)
snp <- snp_attachExtdata()
G <- snp$genotypes
test0 <- parallelVectorSum(G)
G2 <- big_copy(G, ind.row = rep(rows_along(G), 100),
               ind.col = rep(cols_along(G), 10))
dim(G2)

system.time(test1 <- parallelVectorSum(G2))
# 10 sec when parallel over i -> bad
# 1 sec when parallel over j -> very good -> but thread-safe??
str(test1)
# num [1:45420] 48000 13400 16900 12600 34200 24800 9100 33600 36900 31600 ...
*/
