/******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
NumericVector sp_colSumsSq_sym(std::vector<size_t> p,
                               const IntegerVector& i,
                               const NumericVector& x) {

  int m = p.size() - 1;
  NumericVector res(m);

  for (int j = 0; j < m; j++) {

    size_t lo = p[j];
    size_t up = p[j + 1];

    for (size_t k = lo; k < up; k++) {

      int    ind = i[k];
      double val = x[k];

      res[j] += val * val;
      if (ind != j) res[ind] += val * val;  // do not double-count the diagonal
    }
  }

  return res;
}

/******************************************************************************/
