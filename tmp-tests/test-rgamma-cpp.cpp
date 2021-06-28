#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector test_rgamma(int n, double a, double b) {

  NumericVector res(n);
  for (int i = 0; i < n; i++)
    res[i] = ::Rf_rgamma(a, b);

  return res;
}

/*** R
hist(test_rgamma(1e5, 2.3, 5.6))
*/
