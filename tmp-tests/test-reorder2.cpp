/******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

inline double square(double x) { return x * x; }

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
      double r2_ij = x[k];
      score += r2_ij * square(pos[i[k]] - pos_j);
    }
  }

  return score;
}

/******************************************************************************/

double best_switch(const std::vector<size_t>& p,
                   const IntegerVector& i,
                   const NumericVector& x,
                   IntegerVector& pos,
                   std::vector< std::pair<int,double> >& to_consider,
                   double current_best_score,
                   bool sample_best,
                   int max_dist,
                   int j) {

  int best_ind = j;
  double best_delta = 0;
  to_consider.clear();

  size_t lo = p[j];
  size_t up = p[j + 1];
  int pos_j = pos[j];

  // search for a better switch (pos[j] <--> pos[j2])
  int m = pos.size();
  for (int j2 = 0; j2 < m; j2++) {

    if (j2 == j) continue;  // no change

    int pos_j2 = pos[j2];
    if (std::abs(pos_j2 - pos_j) > max_dist) continue;  // too far away

    // vec_diff2 <- (pos[k1] - pos)^2 - (pos[k2] - pos)^2
    // delta3 <- sum(((corr2[, k2] - corr2[, k1]) * vec_diff2)[-ind])
    double delta = 0;
    for (size_t k = lo; k < up; k++) {

      int i_k = i[k];
      if (i_k == j || i_k == j2) continue;

      double r2_ij = x[k];
      double pos_i = pos[i_k];

      delta += r2_ij * (square(pos_j2 - pos_i) - square(pos_j - pos_i));
    }

    size_t lo2 = p[j2], up2 = p[j2 + 1];
    for (size_t k = lo2; k < up2; k++) {

      int i_k = i[k];
      if (i_k == j || i_k == j2) continue;

      double r2_ij2 = x[k];
      double pos_i = pos[i_k];

      delta += r2_ij2 * (square(pos_j - pos_i) - square(pos_j2 - pos_i));
    }

    if (sample_best && delta < 0) {
      to_consider.push_back(std::make_pair(j2, delta));
    }
    if (delta < best_delta) {
      // store the better solution
      best_delta = delta;
      best_ind = j2;
    }
  }

  if (best_ind != j) {

    if (sample_best) {
      int L = to_consider.size();
      NumericVector probs(L);
      for (int l = 0; l < L; l++) probs[l] = -to_consider[l].second;
      int w = Rcpp::sample(L, 1, true, probs, false)[0];
      best_ind   = to_consider[w].first;
      best_delta = to_consider[w].second;
    }

    // use the best solution found -> apply the switch
    pos[j] = pos[best_ind];
    pos[best_ind] = pos_j;
    return current_best_score + 2 * best_delta;

  } else {
    return current_best_score;
  }
}

/******************************************************************************/

// [[Rcpp::export]]
IntegerVector reorder2(const S4& corr,
                       int max_dist,
                       bool sample_index,
                       bool sample_best,
                       int nb_iter = 20) {

  std::vector<size_t> p = corr.slot("p");
  IntegerVector i = corr.slot("i");
  NumericVector x = corr.slot("x");

  int m = p.size() - 1;
  IntegerVector pos = Rcpp::seq_len(m);
  std::vector< std::pair<int,double> > to_consider;

  double score = get_score(p, i, x, pos);
  Rcout << "Initial score: " << score << std::endl;

  for (int iter = 0; iter < nb_iter; iter++) {

    double prev_score = get_score(p, i, x, pos);

    if (sample_index) {
      for (int rep = 0; rep < m; rep++) {
        int j = Rcpp::sample(m, 1, true, R_NilValue, false)[0];
        score = best_switch(p, i, x, pos,
                            to_consider, score, sample_best, max_dist, j);
      }
    } else {
      for (int j = 0; j < m; j++) {
        score = best_switch(p, i, x, pos,
                            to_consider, score, sample_best, max_dist, j);
      }
    }

    double new_score = get_score(p, i, x, pos);
    Rcout << "Iteration #" << iter + 1 << ": " << score <<
      " // " << new_score << std::endl;

    if (new_score == prev_score) break;
  }

  pos.attr("score") = score;
  return pos;
}

/******************************************************************************/
