#include <Rcpp.h>
using namespace Rcpp;

// Based on code from R package {qap}

// [[Rcpp::export]]
IntegerVector qapSA2(const NumericMatrix& A,
                     const NumericMatrix& B,
                     const IntegerVector& perm0,
                     double fiter = 1.1,
                     double ft = 0.5,
                     int maxstep = 50) {

  IntegerVector perm = Rcpp::clone(perm0), best_perm = Rcpp::clone(perm0);

  int n = A.nrow();
  int m = 2 * n;

  NumericVector probs(n);

  // double t = sum(A) * sum(B) / (n * double(n) - n);
  double ol = 0;
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      ol += A(i, j) * B(perm[i], perm[j]);
  double t = (ol / n) / n;
  Rcout << t << std::endl;
  t = 1e-4;

  double min = ol, max = R_NegInf;
  int i1, i2, ibild, jbild, kbild;

  for (int step = 0; step < maxstep; step++) {

    for (int i = 0; i < m; i++) {

      Rcpp::checkUserInterrupt();

      // sample one index
      i1 = n * unif_rand();
      ibild = perm[i1];
      // sample a second one (different from the first one)
      probs = 1.0 / abs(perm - ibild);
      probs[i1] = 0;
      // do { i2 = n * unif_rand(); } while (i2 == i1);
      i2 = sample(n, 1, true, probs, false)[0];
      jbild = perm[i2];

      // evaluate change in obj. function value corresponding to
      // transposition of perm(i1) and perm(i2)
      // NOTE: the following only works for symmetric A and B
      double delta = 0;
      for (int j1 = 0; j1 < n; j1++) {
        if (j1 == i1 || j1 == i2) continue;
        kbild = perm[j1];
        delta -= (A(i1, j1) - A(i2,j1)) * (B(ibild, kbild) - B(jbild, kbild));
      }
      delta = 2 * delta - (A(i1, i1) - A(i2, i2)) * (B(ibild, ibild) - B(jbild, jbild));

      // if obj. function value decreases, always perform transposition
      if (delta > 0) {  // if not, give it a chance to be accepted
        double dt = delta / t;
        if (dt > 10 || unif_rand() > ::exp(-dt)) continue;  // rejected
      }

      // accept transposition and set new permutation perm
      perm[i1] = jbild;
      perm[i2] = ibild;
      ol += delta;

      // adjust bounds for stopping criterion
      if (ol < min) min = ol;
      if (ol > max) max = ol;
      if (max <= min) break;
    }

    // adjust iteration control variables m and t
    m *= fiter;
    t *= ft;

    Rcout << "Iteration #" << step + 1 << ": " << ol << std::endl;
  }

  return perm;
}
