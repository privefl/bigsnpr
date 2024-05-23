/******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

/******************************************************************************/

// Expect both vectors to be sorted
// [[Rcpp::export]]
bool any_overlap(const std::vector<size_t>& P,
                 const std::vector<int>& I,
                 int j1, int j2) {

  size_t lo1 = P[j1], up1 = P[j1 + 1];
  size_t lo2 = P[j2], up2 = P[j2 + 1];
  if (lo1 == up1 || lo2 == up2) return false;

  size_t k1 = lo1, k2 = lo2;
  int i1 = I[k1], i2 = I[k2];

  while (true) {

    if (i1 == i2) return true;      // found one!

    if (i1 < i2) {
      k1++;
      if (k1 == up1) return false;  // no more indices to consider
      i1 = I[k1];
    } else {
      k2++;
      if (k2 == up2) return false;  // no more indices to consider
      i2 = I[k2];
    }
  }
}

/******************************************************************************/

// [[Rcpp::export]]
ListOf<IntegerVector> find_indirect_corr(const std::vector<size_t>& P,
                                         const std::vector<int>& I,
                                         int ncores) {

  int m = P.size() - 1;
  List res(m);
  int chunk_size = ceil(m / (10.0 * ncores));

  #pragma omp parallel num_threads(ncores)
  {
    std::vector<int> keep_with_j0(m);

    #pragma omp for schedule(dynamic, chunk_size)
    for (int j0 = 0; j0 < m; j0++) {

      keep_with_j0.clear();

      for (int j = 0; j < m; j++) {
        if (any_overlap(P, I, j0, j))
          keep_with_j0.push_back(j + 1);
      }

      #pragma omp critical
      res[j0] = wrap(keep_with_j0);
    }
  }

  return res;
}

/******************************************************************************/
