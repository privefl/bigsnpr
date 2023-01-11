/******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
inline int update_pos(int pos, int prev_pos, int new_pos) {

  if (pos == prev_pos) return new_pos;

  if (pos > prev_pos) pos--;
  if (pos >= new_pos) pos++;
  return pos;
}

// [[Rcpp::export]]
double get_new_score(const std::vector<size_t>& p,
                     const IntegerVector& i,
                     const NumericVector& x,
                     const IntegerVector& pos,
                     int prev_pos,
                     int new_pos,
                     double best_score) {

  int m = p.size() - 1;
  IntegerVector pos2(m);
  for (int k = 0; k < m; k++) pos2[k] = update_pos(pos[k], prev_pos, new_pos);

  double score = 0;
  for (int j = 0; j < m; j++) {

    size_t lo = p[j];
    size_t up = p[j + 1];

    int pos2_j = pos2[j];
    for (size_t k = lo; k < up; k++) {
      double dist_ij = pos2[i[k]] - pos2_j;
      double r2_ij = x[k];
      score += r2_ij * dist_ij * dist_ij;
      if (score > best_score) break;  // won't be better anyway, so can stop now
    }
  }

  return score;
}

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

double find_best_place(const std::vector<size_t>& p,
                       const IntegerVector& i,
                       const NumericVector& x,
                       IntegerVector& pos,
                       double best_score,
                       int max_dist,
                       int j) {

  int pos_j = pos[j];
  int best_pos = pos_j;

  // try to find a better position for j
  int m = pos.size();
  for (int new_pos_j = 0; new_pos_j < m; new_pos_j++) {

    if (new_pos_j == pos_j) continue;  // no change
    if (std::abs(new_pos_j - pos_j) > max_dist) continue;  // too far apart

    double score = get_new_score(p, i, x, pos, pos_j, new_pos_j, best_score);

    // found a better one
    if (score < best_score) {
      best_score = score;
      best_pos = new_pos_j;
    }
  }

  // found a better one -> apply the change
  if (best_pos != pos_j) {
    for (int k = 0; k < m; k++) pos[k] = update_pos(pos[k], pos_j, best_pos);
  }

  return best_score;
}

/******************************************************************************/

// [[Rcpp::export]]
IntegerVector reorder3(const S4& corr, int max_dist, int nb_iter = 20) {

  std::vector<size_t> p = corr.slot("p");
  IntegerVector i = corr.slot("i");
  NumericVector x = corr.slot("x");

  int m = p.size() - 1;
  IntegerVector pos = Rcpp::seq_len(m);

  double score = get_score(p, i, x, pos);
  Rcout << "Initial score: " << score << std::endl;

  for (int iter = 0; iter < nb_iter; iter++) {

    double prev_score = score;

    for (int j = 0; j < m; j++) {
      score = find_best_place(p, i, x, pos, score, max_dist, j);
    }

    Rcout << "Iteration #" << iter + 1 << ": " << score <<
      " // " << get_score(p, i, x, pos) << std::endl;

    if (score > (0.999 * prev_score)) break;
  }

  return pos;
}

/******************************************************************************/
