/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
List get_L(const std::vector<size_t>& p,
           const IntegerVector& i,
           const NumericVector& x,
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

double square(double x) { return x * x; }

// [[Rcpp::export]]
List get_C(const arma::sp_mat& L, const NumericVector& pos,
           double min_size, double max_size, int K) {

  int m = L.n_rows;  // L now has an extra column with all 0s for convenience
  std::vector< std::vector<float> > res_E(m);

  IntegerVector first_col(m, NA_INTEGER);

  for (int col = 0; col < m; col++) {

    double e = 0;
    double min_pos = pos[col] - max_size;
    double max_pos = pos[col] - min_size;

    for (int row = col; row >= 0; row--) {

      e += L(row, col + 1);

      double pos_row = pos[row];

      if (pos_row < min_pos) break;
      if (pos_row <= max_pos) {
        // E[row, col] = e (can deduce 'col' from 'row' if store first)
        if (res_E[row].size() == 0) first_col[row] = col;
        res_E[row].push_back(e);
      }
    }
  }

  // Computing all minimum costs and corresponding indices
  IntegerMatrix best_ind(m, K); best_ind.fill(NA_INTEGER);
  NumericMatrix C(m, K); C.fill(NA_REAL);

  // Only a few indices allow one block only
  double min_pos = pos[m - 1] - max_size;
  double max_pos = pos[m - 1] - min_size;
  for (int row = m - 1; row >= 0; row--) {
    double pos_row = pos[row];
    if (pos_row < min_pos) break;
    if (pos_row <= max_pos) {
      best_ind(row, 0) = m;
      C(row, 0) = 0;
    }
  }

  // Iterating over total numbers of blocks allowed
  for (int k = 1; k < K; k++) {

    // starting from the end
    for (int row = m - 2; row >= 0; row--) {

      int col = first_col[row];
      std::vector<float> E_i = res_E[row];
      double min_cost = R_PosInf;

      for (auto it = E_i.begin(); it != E_i.end(); ++it) {
        double cost = square(*it) + C(col + 1, k - 1);
        if (cost < min_cost) {
          min_cost = cost;
          best_ind(row, k) = col + 1;
          C(row, k) = cost;
        }
        col++;
      }
    }
  }

  return List::create(_["C"] = C, _["best_ind"] = best_ind);
}

/******************************************************************************/

// // [[Rcpp::export]]
// double test_pow(int n) {
//   return std::pow(n, 1.5);
// }

// // [[Rcpp::export]]
// IntegerVector test_seq() {
//   return seq(1, 1);
// }

// // [[Rcpp::export]]
// bool test_compare_NA() {
//   double cost = NA_REAL;
//   return (cost < R_PosInf);
// }
