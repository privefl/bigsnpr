#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double dbeta_cpp(double x, double a, double b) {
  return ::Rf_dbeta(x, a, b, 0);
}
// [[Rcpp::export]]
double dgamma_cpp(double x, double a, double b) {
  return ::Rf_dgamma(x, a, b, 0);
}


/*** R
a <- runif(1)
b <- runif(1)
x <- runif(1)
dbeta_cpp(x, a, b)
dbeta(x, a, b)
dgamma_cpp(x, a, b)
dgamma(x, a, b)
dgamma(x, a, 1 / b)
dgamma(x, a, scale = b)
*/
