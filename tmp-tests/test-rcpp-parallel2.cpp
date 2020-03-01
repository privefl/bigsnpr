// [[Rcpp::depends(RcppParallel, rmio, bigstatsr)]]
#include <bigstatsr/BMCodeAcc.h>
#include <RcppParallel.h>
using namespace RcppParallel;

struct Sum : public Worker {

  SubBMCode256Acc macc;
  std::size_t j0, j;

  // accumulated xySum
  double xySum;

  // constructors
  Sum(SubBMCode256Acc macc, std::size_t j0, std::size_t j) :
    macc(macc), j0(j0), j(j), xySum(0) {}
  Sum(const Sum& sum, Split) : macc(sum.macc), j0(sum.j0), j(sum.j), xySum(0) {}

  // accumulate just the element of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      xySum += macc(i, j) * macc(i, j0);
    }
  }

  // join my xySum with that of another Sum
  void join(const Sum& rhs) {
    xySum += rhs.xySum;
  }
};

// [[Rcpp::export]]
NumericVector parallelVectorSum(Environment BM) {

  XPtr<FBM> xpBM = BM["address"];
  std::size_t n = xpBM->nrow();
  std::size_t m = xpBM->ncol();
  SubBMCode256Acc macc(xpBM, seq_len(n), seq_len(m), BM["code256"], 1);

  NumericVector res(m);
  for (size_t j = 0; j < m; j++) {
    Sum sum(macc, 0, j);
    parallelReduce(0, n, sum);
    res[j] = sum.xySum;
  }

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
system.time(test1 <- parallelVectorSum(G2)) # 42 sec -> very bad
*/
