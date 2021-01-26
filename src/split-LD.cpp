/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
IntegerVector min_row(std::vector<size_t> p, IntegerVector i) {

  int m = p.size() - 1;
  IntegerVector min_row(m);

  for (int col = 0; col < m; col++) min_row[col] = i[p[col]];

  return min_row;
}

/******************************************************************************/

// [[Rcpp::export]]
DataFrame get_L(std::vector<size_t> p, IntegerVector i, NumericVector x,
                double thr_r2) {

  std::vector<int>    res_i;
  std::vector<int>    res_j;
  std::vector<double> res_x;

  int m = p.size() - 1;

  for (int col = 0; col < m; col++) {

    double l = 0;
    size_t k = p[col + 1] - 1; // safe because corr(col, col) = 1

    for (int row = i[k]; row > col; row--) {

      if (row == i[k]) {
        double r2 = x[k] * x[k];
        if (r2 >= thr_r2) l += r2;
        k--;
      }

      if (l > 0) {
        res_i.push_back(col);
        res_j.push_back(row);
        res_x.push_back(l);
      }
    }
  }

  return DataFrame::create(_["i"] = res_i, _["j"] = res_j, _["x"] = res_x);
}

/******************************************************************************/

// [[Rcpp::export]]
DataFrame get_E(const arma::sp_mat& L, IntegerVector min_row) {

  std::vector<int> res_i;
  std::vector<int> res_j;
  std::vector<double> res_x;

  int m2 = min_row.size() - 1;

  for (int col = 0; col < m2; col++) {

    int last = min_row[col];
    double e = 0;

    for (int row = col; row >= last; row--) {
      e += L(row, col + 1);
      res_i.push_back(row + 1);
      res_j.push_back(col + 1);
      res_x.push_back(e);
    }
  }

  return DataFrame::create(_["i"] = res_i, _["j"] = res_j, _["x"] = res_x);
}

/******************************************************************************/
