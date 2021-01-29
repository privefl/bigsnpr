/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
List get_L(std::vector<size_t> p, IntegerVector i, NumericVector x,
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

  return List::create(_["i"] = res_i, _["j"] = res_j, _["x"] = res_x);
}

/******************************************************************************/

// [[Rcpp::export]]
List get_C(const arma::sp_mat& L, int min_size, int max_size, double lambda) {

  int m = L.n_rows;  // L now has an extra column with all 0s for convenience
  std::vector< std::vector<float> > res_E(m);

  for (int col = 0; col < m; col++) {

    double e = 0;
    int count = 0;

    for (int row = col; row >= 0; row--) {

      e += L(row, col + 1);
      count++;

      if (count >= min_size) {
        // E[row, col] = e with reparametrization to save memory
        // (can deduce 'col' from 'row' and 'min_size')
        res_E[row].push_back(e);
        if (count == max_size) break;
      }
    }
  }

  // Precomputing penalization costs
  NumericVector pena_cost(m);
  for (int diff = 0; diff < m; diff++)
    pena_cost[diff] = lambda * std::max(0, diff + 1 - min_size);

  // Computing all minimum costs, starting from the end
  IntegerVector best_ind(m, NA_INTEGER);
  NumericVector C(m + 1, NA_REAL); C[m] = 0;
  // cannot split -> need to be used as a single block
  int last_split_once = std::min(2 * min_size - 1, max_size);
  C[m - seq(min_size, last_split_once)] = 0;
  best_ind[m - seq(min_size, last_split_once)] = m;
  // for (int row = m - max_size - 1; row >= 0; row--) {
  for (int row = m - last_split_once - 1; row >= 0; row--) {
    std::vector<float> E_i = res_E[row];
    // Rcout << E_i.size() << std::endl;
    double min_cost = R_PosInf;
    int col = row + min_size - 1;
    for (auto it = E_i.begin(); it != E_i.end(); ++it) {
      double cost = double(*it) + C[col + 1] + pena_cost[col - row];
      if (cost < min_cost) {
        best_ind[row] = col + 1;
        C[row] = cost;
        min_cost = cost;
      }
      col++;
    }
  }

  return List::create(_["C"] = C, _["best_ind"] = best_ind);
}

/******************************************************************************/

// // [[Rcpp::export]]
// double test_pow(int n) {
//   return std::pow(n, 1.5);
// }
//
// // [[Rcpp::export]]
// IntegerVector test_seq() {
//   return seq(1, 1);
// }
