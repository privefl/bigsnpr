/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>

using namespace Rcpp;

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
List get_C(const arma::sp_mat& L, int min_size, int max_size, double lambda) {

  int m = L.n_cols;
  int m2 = m - 1;

  std::vector< std::vector<float> > res_E(m);

  for (int col = 0; col < m2; col++) {

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

  // Computing all minimum costs, starting from the end
  IntegerVector best_ind(m, NA_INTEGER);
  NumericVector C(m, NA_REAL);
  C[m - seq(min_size, max_size)] = 0;
  for (int row = m - max_size - 1; row >= 0; row--) {
    std::vector<float> E_i = res_E[row];
    double min_cost = R_PosInf;
    int col = row + min_size - 1;
    for (auto it = E_i.begin(); it != E_i.end(); ++it) {
      double cost = *it + C[col + 1] + lambda * (col - row);
      if (cost < min_cost) {
        best_ind[row] = col;
        C[row] = *it + C[col + 1];
        min_cost = cost;
      }
      col++;
    }
  }

  return List::create(_["C"] = C, _["best_ind"] = best_ind + 1);
}

/******************************************************************************/
