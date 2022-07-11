/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>

using namespace Rcpp;

/******************************************************************************/

inline double square(double x) {
  return x * x;
}

/******************************************************************************/

// [[Rcpp::export]]
List get_L(std::vector<size_t> p,
           const IntegerVector& i,
           const NumericVector& x,
           double thr_r2,
           double max_r2) {

  // L(i, j) = sum_{q=j}^m r(i, q)^2 (note: start at j, but j+1 in the paper)
  // L(i, i) is not used in E(i, j), so not computed

  std::vector<int>    res_i;
  std::vector<int>    res_j;
  std::vector<double> res_x;

  int m = p.size() - 1;

  for (int col = 0; col < m; col++) {

    double l = 0;
    size_t k = p[col + 1] - 1; // safe because corr(col, col) = 1

    // start from the last non-zero element of the column
    // stop at col+1
    for (int row = i[k]; row > col; row--) {

      if (row == i[k]) {
        double r2 = x[k] * x[k];
        if (r2 >= thr_r2) {
          if (r2 > max_r2) {
            l = R_PosInf;
          } else {
            l += r2;
          }
        }
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
List get_C(const arma::sp_mat& L, int min_size, int max_size, int max_K, double max_cost) {

  int m = L.n_rows;  // L now has an extra column with all 0s for convenience
  std::vector< std::vector<float> > res_E(m);

  // E(i, j) = sum_{p=i}^j L(p, j+1)

  for (int col = 0; col < m; col++) {

    double e = 0;
    int count = 0;

    for (int row = col; row >= 0; row--) {

      e += L(row, col + 1);  // compute E(j, j), then E(j-1, j), etc
      if (e > max_cost) break;
      count++;

      if (count >= min_size) {
        // E(row, col) = e with reparametrization to save memory
        // -> the first row has to be (col - min_size + 1)
        res_E[col].push_back(e);
        if (count == max_size) break;
      }
    }
  }

  // Computing all minimum costs and corresponding indices
  IntegerMatrix best_ind(m, max_K); best_ind.fill(NA_INTEGER);
  NumericMatrix C1(m, max_K); C1.fill(R_PosInf);
  NumericMatrix C2(m, max_K); C2.fill(R_PosInf);

  // Only a few indices allow one block only
  for (auto size : seq(min_size, max_size)) {
    best_ind(m - size, 0) = m;
    C1(m - size, 0) = 0;
    C2(m - size, 0) = square(size);
  }

  // Iterating over total numbers of blocks allowed
  for (int k = 1; k < max_K; k++) {

    // starting from the end
    for (int col = m - 1; col >= 0; col--) {

      std::vector<float> E_j = res_E[col];
      int row = col - min_size + 1;

      for (auto it = E_j.begin(); it != E_j.end(); it++, row--) {
        double cost1 = double(*it) + C1(col + 1, k - 1);
        if (cost1 < C1(row, k)) {
          best_ind(row, k) = col + 1;
          C1(row, k) = cost1;
          C2(row, k) = square(col - row + 1) + C2(col + 1, k - 1);
        } else if (cost1 == C1(row, k)) {
          double cost2 = square(col - row + 1) + C2(col + 1, k - 1);
          if (cost2 < C2(row, k)) {
            best_ind(row, k) = col + 1;
            C2(row, k) = cost2;
          }
        }
      }
    }

    if (C1(0, k) > max_cost && C1(0, k) > C1(0, k - 1)) break;
  }

  return List::create(_["C"] = C1, _["best_ind"] = best_ind);
}

/******************************************************************************/

// [[Rcpp::export]]
double get_perc(const NumericVector& p,  // integers, but possibly very large
                const IntegerVector& i,
                const IntegerVector& all_last) {

  int m = p.size() - 1;
  double count_all = 2 * i.size() - m;  // count the diaginal only once
  double count_within = count_all;

  int grp_num = 0;
  int limit = all_last[grp_num];

  for (int j = 0; j < m; j++) {

    if (j > limit) {
      grp_num++;
      limit = all_last[grp_num];
    }

    size_t lo = p[j];      // the diagonal -> can skip it
    size_t up = p[j + 1];  // up > lo because of the diag

    // start by the end to stop early
    for (size_t k = up - 1; k > lo; k--) {
      if (i[k] > limit) {
        count_within -= 2;
      } else {
        break;  // all remaining i are in the block
      }
    }
  }

  return count_within / count_all;
}

/******************************************************************************/
