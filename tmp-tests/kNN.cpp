#include <Rcpp.h>
using namespace Rcpp;


// x should be ordered
// [[Rcpp::export]]
NumericMatrix dist_kNN(const NumericVector& x,
                       const NumericVector& y,
                       int nfirst) {

  int n = x.size();

  NumericVector min_dist(nfirst);
  NumericMatrix res(n, nfirst);

  for (int i = 0; i < n; i++) {

    // Rcout << i << ":";

    for (int k = 0; k < nfirst; k++) min_dist[k] = R_PosInf;
    int pos_max = 0;
    double d_max = R_PosInf;

    bool go_up = true, go_down = true;
    int j, l = 1;

    while (go_up || go_down) {

      j = i + l;
      // Rcout << " " << j;
      go_up &= (j < n);
      if (go_up) {
        double d1 = x[j] - x[i];
        if (d1 > d_max) {
          go_up = false;
        } else {
          double d2 = y[j] - y[i];
          double d = ::sqrt(d1 * d1 + d2 * d2);
          if (d < d_max) { // new max amongst min distances
            min_dist[pos_max] = d;
            pos_max = which_max(min_dist);
            d_max = min_dist[pos_max];
          }
        }
      }

      j = i - l;
      go_down &= (j >= 0);
      if (go_down) {
        double d1 = x[i] - x[j];
        if (d1 > d_max) {
          go_down = false;
        } else {
          double d2 = y[i] - y[j];
          double d = ::sqrt(d1 * d1 + d2 * d2);
          if (d < d_max) { // new max amongst min distances
            min_dist[pos_max] = d;
            pos_max = which_max(min_dist);
            d_max = min_dist[pos_max];
          }
        }
      }

      l++;
    }

    // Rcout << std::endl;

    std::sort(min_dist.begin(), min_dist.end());
    for (int k = 0; k < nfirst; k++) res(i, k) = min_dist[k];
  }

  return res;
}
