/******************************************************************************/

#include <Rcpp.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
List get_L(std::vector<size_t> p, IntegerVector i, NumericVector x) {

  std::vector<int> res_i;
  std::vector<int> res_j;
  std::vector<double> res_x;

  int m = p.size() - 1;
  Rcout << m << std::endl;

  for (int j = 0; j < m; j++) {

    size_t lo = p[j];
    size_t up = p[j + 1];

    if (up > lo) {
      double l = 0;
      for (size_t k = up - 1; ;k--) {
        // Rcout << k << std::endl;
        double val = x[k];
        l += val * val;
        res_i.push_back(j);
        res_j.push_back(i[k]);
        res_x.push_back(l);
        if (k == lo) break;
      }
    }
  }

  return List::create(res_i, res_j, res_x);
}

/******************************************************************************/

// [[Rcpp::export]]
List get_E(std::vector<size_t> p, IntegerVector i, NumericVector x) {

  std::vector<int> res_i;
  std::vector<int> res_j;
  std::vector<double> res_x;

  int m2 = p.size() - 2;
  Rcout << m2 << std::endl;

  for (int j = 0; j < m2; j++) {

    size_t lo = p[j + 1];
    size_t up = p[j + 2];

    if (up > lo) {
      double l = 0;
      for (size_t k = up - 1; ;k--) {
        // Rcout << k << std::endl;
        l += x[k];
        res_i.push_back(i[k] + 1);
        res_j.push_back(j + 1);
        res_x.push_back(l);
        if (k == lo) break;
      }
    }
  }

  return DataFrame::create(_["i"] = res_i, _["j"] = res_j, _["x"] = res_x);
}

/******************************************************************************/

// [[Rcpp::export]]
List get_E2(std::vector<size_t> p, IntegerVector i, NumericVector x) {

  std::vector<int> res_i;
  std::vector<int> res_j;
  std::vector<double> res_x;

  int m2 = p.size() - 2;
  Rcout << m2 << std::endl;

  ListOf< std::vector<int> > res(m2);

  for (int j = 0; j < m2; j++) {

    size_t lo = p[j + 1];
    size_t up = p[j + 2];

    if (up > lo) {
      double l = 0;
      for (size_t k = up - 1; ;k--) {
        // Rcout << k << std::endl;
        l += x[k];
        res_i.push_back(i[k] + 1);
        res_j.push_back(j + 1);
        res_x.push_back(l);
        if (k == lo) break;
      }
    }
  }

  return res;
}

/******************************************************************************/
