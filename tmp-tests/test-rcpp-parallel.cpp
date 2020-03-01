// [[Rcpp::depends(RcppParallel, rmio, bigstatsr)]]
#include <bigstatsr/BMCodeAcc.h>
#include <RcppParallel.h>
using namespace RcppParallel;

struct Sum : public Worker {

  SubBMCode256Acc macc;
  RVector<double> output;
  std::size_t j0, n;

  // constructors
  Sum(SubBMCode256Acc macc, std::size_t j0, NumericVector output) :
    macc(macc), output(output), j0(j0), n(macc.nrow()) {}

  // accumulate just the element of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++) {
      double xySum = 0;
      for (std::size_t i = 0; i < n; i++) {
        xySum += macc(i, j) * macc(i, j0);
      }
      output[j] = xySum;
    }
  }
};

// [[Rcpp::export]]
NumericVector parallelVectorSum(Environment BM) {

  XPtr<FBM> xpBM = BM["address"];
  std::size_t n = xpBM->nrow();
  std::size_t m = xpBM->ncol();
  SubBMCode256Acc macc(xpBM, seq_len(n), seq_len(m), BM["code256"], 1);

  NumericVector res(m);
  Sum sum(macc, 0, res);
  parallelFor(0, m, sum, 1);

  return res;
}


/*** R
RcppParallel::setThreadOptions(2)
library(bigsnpr)
snp <- snp_attachExtdata()
G <- snp$genotypes
test0 <- parallelVectorSum(G)
G2 <- big_copy(G, ind.row = rep(rows_along(G), 100),
               ind.col = rep(cols_along(G), 10))
dim(G2)

RcppParallel::setThreadOptions(1)
system.time(test1 <- parallelVectorSum(G2))
system.time(test1 <- parallelVectorSum(G2)) # 5.4 / 5.0
RcppParallel::setThreadOptions(2)
system.time(test2 <- parallelVectorSum(G2)) # 3.4 / 2.1
RcppParallel::setThreadOptions(3)
system.time(parallelVectorSum(G2))          # 1.6 / 1.6
RcppParallel::setThreadOptions(4)
system.time(parallelVectorSum(G2))          # 1.4 / 1.3
*/
