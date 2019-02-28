#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector sample_from_prob(double p0, double p1, int n) {
  IntegerVector tab(3);
  for (int i = 0; i < n; i++) {
    double first = ::unif_rand() - p0;
    int g = (first < 0) ? 0 : ((first < p1) ? 1 : 2);
    (tab[g])++;
  }
  return tab;
}

/*** R
p0 <- 0.1
p1 <- 0.9

system.time({
  gen <- runif(1e7) - p0
  print(table(ifelse(gen < 0, 0L, ifelse(gen < p1, 1L, 2L))))
})

system.time(print(sample_from_prob(0.5, 0.2, 1e9)))
*/
