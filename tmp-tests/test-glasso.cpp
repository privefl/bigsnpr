#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double abs_dist(const NumericVector& X,
                const NumericVector& Y) {

  int K = X.size();
  double dist = 0;
  for (int k = 0; k < K; k++)
    dist += std::abs(X[k] - Y[k]);

  return dist;
}

// [[Rcpp::export]]
void update_WXj(NumericVector& WXj,
                const NumericMatrix& W,
                int i,
                double delta) {
  // WXj = WXj + W[, i] * delta
  int K = WXj.size();
  for (int k = 0; k < K; k++)
    WXj[k] += W(k, i) * delta;
}

inline double sgn(double x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

// [[Rcpp::export]]
double inner_loop(const NumericMatrix& S,
                  NumericMatrix& W,
                  NumericMatrix& X,
                  NumericVector& WXj,
                  const NumericVector& Wd,
                  int n,
                  double rho,
                  double thrLasso) {

  // int count = 0;

  double dw = 0;
  for (int j = 0; j < n; j++) {

    WXj.fill(0);
    for (int i = 0; i < n; i++) {
      if (X(i, j) != 0) {
        for (int k = 0; k < n; k++)
          WXj[k] = WXj[k] + W(k, i) * X(i, j);
      }
    }

    while (true) {

      double dlx = 0;
      for (int i = 0; i < n; i++) {
        if (i != j) {
          double a = S(i, j) - WXj[i] + Wd[i] * X(i, j);
          double b = std::abs(a) - rho;
          double c = (b > 0) ? sgn(a) * b / Wd[i] : 0;
          double delta = c - X(i, j);
          if (delta != 0) {
            X(i, j) = c;
            dlx = std::max(dlx, std::abs(delta));
            // count++;
            for (int k = 0; k < n; k++)
              WXj[k] += W(k, i) * delta;
          }
        }
      }

      if (dlx < thrLasso) break;
    }

    WXj[j] = Wd[j];
    double dist = 0;
    for (int k = 0; k < n; k++)
      dist += std::abs(WXj[k] - W(k, j));
    dw = std::max(dw, dist);
    for (int k = 0; k < n; k++) {
      // if (S(k, j) != 0) {
      //   W(k, j) = WXj[k];
      //   W(j, k) = W(k, j);
      // }
      W(k, j) = WXj[k];
      W(j, k) = W(k, j);
    }
  }

  return dw;
}



// [[Rcpp::export]]
NumericVector test_abs(const NumericVector& X) {
  return abs(X);
}

// [[Rcpp::export]]
double test_abs2(double x) {
  return abs(x);
}

