/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <bigsnpr/bed-acc.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
void read_plink(Environment BM,
                const std::string filename, int n_ind, int n_snp,
                const RawMatrix& decode) {

  // Init FBM accessor
  XPtr<FBM> xpBM = BM["address"];
  BMAcc<unsigned char> macc(xpBM);

  // Init bed accessor
  bedAcc bacc(filename, n_ind, n_snp, decode);

  size_t i, j, n = n_ind, m = n_snp;
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      macc(i, j) = bacc(i, j);
    }
  }
}

/******************************************************************************/
