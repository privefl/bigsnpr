#include <Rcpp.h>
using namespace Rcpp;
#include <bigsnpr/SubTripleAcc.h>

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
ListOf<NumericVector> timesTwo(NumericVector x) {
  return List::create(x, x > 0.5, x > 1.5);
}

// [[Rcpp::export]]
int test(int x) {
  return x / 3;
}



template <class C>
ListOf<NumericVector> bigcolvars(C macc) {
  int n = macc.nrow();
  int m = macc.ncol();

  NumericVector res(m), res2(m);
  double x, xSum, xxSum;
  int i, j;

  for (j = 0; j < m; j++) {
    xSum = xxSum = 0;
    for (i = 0; i < n; i++) {
      x = macc(i, j);
      xSum += x;
      xxSum += x*x;
    }
    res[j] = xxSum - xSum * xSum / n;
    res2[j] = xSum;
  }

  return List::create(_["sum"] = res2,
                      _["var"] = res/(n-1));
}

// Dispatch function for bigcolvars
// [[Rcpp::export]]
ListOf<NumericVector> bigcolvars(const S4& BM,
                                 const IntegerVector& rowInd,
                                 const IntegerVector& colInd) {

  XPtr<BigMatrix> xpMat = BM.slot("address");
  IntegerVector rows = rowInd - 1;
  IntegerVector cols = colInd - 1;

  //Rf_inherits(BM, "BM.code"))
  return bigcolvars(RawSubTripleAcc(*xpMat, rows, cols, BM.slot("code")));
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
(x <- runif(10, 0, 2))
timesTwo(x)

test(5)

desc <- snp_attachExtdata()$genotypes
G <- attach.BM(desc)
test <- bigcolvars(G, rows_along(G), cols_along(G))
true <- big_colstats(G)
*/
