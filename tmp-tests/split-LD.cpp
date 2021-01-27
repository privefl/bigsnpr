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

// // [[Rcpp::export]]
// DataFrame get_E(const arma::sp_mat& L, IntegerVector min_row, int min_size) {
//
//   std::vector<int> res_i;
//   std::vector<int> res_j;
//   std::vector<double> res_x;
//
//   int m2 = min_row.size() - 1;
//
//   for (int col = 0; col < m2; col++) {
//
//     int last = min_row[col];
//     double e = 0;
//     int count = 0;
//
//     for (int row = col; row >= last; row--) {
//       e += L(row, col + 1);
//       if (++count >= min_size) {
//         res_i.push_back(row + 1);
//         res_j.push_back(col + 1);
//         res_x.push_back(e);
//       }
//     }
//   }
//
//   return DataFrame::create(_["i"] = res_i, _["j"] = res_j, _["x"] = res_x);
// }

/******************************************************************************/

// [[Rcpp::export]]
List get_E2(const arma::sp_mat& L, IntegerVector min_row, int min_size) {

  int m = min_row.size();
  int m2 = m - 1;

  std::vector< std::pair< std::vector<int>, std::vector<double> > > res(m);

  // if can write this as a matrix, then can compute in parallel
  for (int col = 0; col < m2; col++) {

    int last = min_row[col];
    double e = 0;
    int count = 0;

    for (int row = col; row >= last; row--) {
      e += L(row, col + 1);
      if (++count >= min_size) {
        res[row].first .push_back(col + 1);
        res[row].second.push_back(e);
      }
    }
  }

  // transform to some R structure
  List res2(m);
  for (int row = m2; row >= 0; row--) {
    res2[row] = List::create(_["j"] = res[row].first,
                             _["x"] = res[row].second);
    res.pop_back(); // free
  }

  return res2;
}

/******************************************************************************/

// [[Rcpp::export]]
List get_E3(const arma::sp_mat& L, int min_size, int max_size, double lambda) {

  int m = L.n_cols;
  int m2 = m - 1;

  std::vector< std::vector<float> > res_E(m);

  // if can write this as a matrix, then can compute in parallel
  for (int col = 0; col < m2; col++) {

    double e = 0;
    int count = 0;

    for (int row = col; row >= 0; row--) {

      e += L(row, col + 1);
      count++;

      if (count >= min_size) {
        // E[row, col] = e
        res_E[row].push_back(e);
        if (count == max_size) break;
      }
    }
  }

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
