/******************************************************************************/

#include <bigstatsr/BMAcc.h>

using namespace Rcpp;
using std::size_t;

/******************************************************************************/

// [[Rcpp::export]]
void replaceSNP(Environment BM,
                Environment BM2,
                const IntegerVector& rowInd,
                const IntegerVector& colInd) {

  XPtr<FBM_RW> xpBM = BM["address_rw"];
  BMAcc_RW<unsigned char> macc(xpBM);

  XPtr<FBM> xpBM2 = BM2["address"];
  SubBMAcc<unsigned char> macc2(xpBM2, rowInd, colInd, 1);

  myassert_size(macc.nrow(), macc2.nrow());
  myassert_size(macc.ncol(), macc2.ncol());

  size_t n = macc.nrow();
  size_t m = macc.ncol();

  // #pragma omp parallel for num_threads(ncores) -> does NOT improve
  for (size_t j = 0; j < m; j++)
    for (size_t i = 0; i < n; i++)
      macc(i, j) = macc2(i, j);
}

/******************************************************************************/
