#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector replaceNA(NumericVector x) {

  x[is_na(x)] = 3;

  return x;
}

/*** R
x <- c(1, 2, NA, 3.5, rep(NA, 4))
replaceNA(x)
*/
