// [[Rcpp::depends(BH, pcadapt, bigstatsr)]]
#include <pcadapt/bed-acc.h>
#include <bigstatsr/BMAcc.h>

// [[Rcpp::export]]
void fillMat(SEXP xptr_bed, Environment BM) {

  XPtr<FBM> xpBM = BM["address"];
  size_t n = xpBM->nrow();
  size_t m = xpBM->ncol();
  SubBMAcc<unsigned char> macc(xpBM, seq_len(n) - 1, seq_len(m) - 1);

  XPtr<bed> xpBed(xptr_bed);
  bedAcc bedacc(xpBed, seq_len(m));

  for (size_t j = 0; j < m; j++) {
    for (size_t i = 0; i < n; i++) {
      macc(i, j) = bedacc(i, j);
    }
  }
}
