/******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
double get_score(const std::vector<size_t>& p,
                 const IntegerVector& i,
                 const NumericVector& x,
                 const IntegerVector& pos) {

  double score = 0;
  int m = p.size() - 1;

  for (int j = 0; j < m; j++) {

    size_t lo = p[j];
    size_t up = p[j + 1];

    int pos_j = pos[j];
    for (size_t k = lo; k < up; k++) {
      double dist_ij = pos[i[k]] - pos_j;
      double r2_ij = x[k];
      score += r2_ij * dist_ij * dist_ij;
    }
  }

  return score;
}

/******************************************************************************/

// [[Rcpp::export]]
bool which_best_place(const std::vector<size_t>& p,
                      const IntegerVector& i,
                      const NumericVector& x,
                      IntegerVector& pos,
                      int j) {

  size_t lo = p[j];
  size_t up = p[j + 1];

  // get current position as current best
  int pos_j = pos[j];
  int best_pos = pos_j;
  double best_score = 0;
  for (size_t k = lo; k < up; k++) {
    double dist_ij = pos[i[k]] - pos_j;
    double r2_ij = x[k];
    best_score += r2_ij * dist_ij * dist_ij;
  }

  // try to find a better one
  int m = pos.size();
  for (int new_pos_j = 0; new_pos_j < m; new_pos_j++) {

    if (new_pos_j == pos_j) continue;  // no change

    double score = 0;

    for (size_t k = lo; k < up; k++) {

      int pos_i = pos[i[k]];
      if (pos_i == pos_j) continue;  // dist is 0

      pos_i -= (pos_i > pos_j);
      pos_i += (pos_i >= new_pos_j);
      double r2_ij = x[k];
      double dist_ij = pos_i - new_pos_j;

      score += r2_ij * dist_ij * dist_ij;
      if (score > best_score) break;
    }

    // found a better one
    if (score < best_score) {
      best_score = score;
      best_pos = new_pos_j;
    }
  }

  // found a better one
  bool has_changed = (best_pos != pos_j);
  if (has_changed) {
    // update all pos
    for (int i = 0; i < m; i++) {
      pos[i] -= (pos[i] > pos_j);
      pos[i] += (pos[i] >= best_pos);
    }
    pos[j] = best_pos;
  }

  return has_changed;
}

/******************************************************************************/

// [[Rcpp::export]]
IntegerVector reorder(const S4& corr, int nb_iter = 20) {

  std::vector<size_t> p = corr.slot("p");
  IntegerVector i = corr.slot("i");
  NumericVector x = corr.slot("x");

  int m = p.size() - 1;
  IntegerVector pos = Rcpp::seq_len(m);

  Rcout << "Initial score: " << get_score(p, i, x, pos) << std::endl;

  for (int iter = 0; iter < nb_iter; iter++) {

    IntegerVector prev_pos = Rcpp::clone(pos);
    double prev_score = get_score(p, i, x, pos);

    for (int j = 0; j < m; j++) {
      bool has_changed = which_best_place(p, i, x, pos, j);
      // if (has_changed) Rcout << "Moved " << j + 1 << std::endl;
    }

    double score = get_score(p, i, x, pos);
    Rcout << "Iteration #" << iter + 1 << ": " << score << std::endl;
    if (prev_score < score) return prev_pos;
  }

  return pos;
}

/******************************************************************************/
