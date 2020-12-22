#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double test1(NumericVector x) {

  double sum = 0;
  int n = x.size();
  for (int i = 0; i < n; i++) {
    sum += x[i];
  }

  return sum;
}

// [[Rcpp::export]]
double test2(NumericVector x, IntegerVector ind) {

  double sum = 0;
  for (const int& i : ind) {
    sum += x[i];
  }

  return sum;
}
