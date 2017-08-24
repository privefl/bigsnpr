/******************************************************************************/

#include <bigstatsr/BMAcc.h>
#include <bigstatsr/utils.h>
#include <Rcpp.h>

using namespace Rcpp;
using std::size_t;

/******************************************************************************/

// [[Rcpp::export]]
void replaceSNP(Environment BM,
                Environment BM2,
                const IntegerVector& rowInd,
                const IntegerVector& colInd) {

  XPtr<FBM> xpBM = BM["address"];
  BMAcc<unsigned char> macc(xpBM);

  XPtr<FBM> xpBM2 = BM2["address"];
  SubBMAcc<unsigned char> macc2(xpBM2, rowInd - 1, colInd - 1);

  myassert(macc.nrow() == macc2.nrow(), ERROR_DIM);
  myassert(macc.ncol() == macc2.ncol(), ERROR_DIM);

  for (size_t j = 0; j < macc.ncol(); j++)
    for (size_t i = 0; i < macc.nrow(); i++)
      macc(i, j) = macc2(i, j);
}

/******************************************************************************/
