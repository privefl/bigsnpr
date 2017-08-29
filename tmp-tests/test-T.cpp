#include <Rcpp.h>
using namespace Rcpp;

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
double testT(double r, int N = 517) {
  double t = r * sqrt((N - 2) / (1 - r*r));
  return t;
}

// [[Rcpp::export]]
double getNA() {
  return NA_REAL;
}

// [[Rcpp::export]]
bool isna(double x) {
  return R_IsNA(x);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
N <- 517
r <- sqrt(0.2)
t <- r * sqrt((N-2)/(1-r^2))
all.equal(testT(r, N), t)

isna(NA)
isna(NA_real_)
getNA()
getNA() == NA_real_
isna(getNA())
*/
