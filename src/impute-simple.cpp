/******************************************************************************/

#include <bigstatsr/BMAcc.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
void impute(Environment BM,
            const IntegerVector& rowInd,
            const IntegerVector& colInd,
            int method) {

  XPtr<FBM_RW> xpBM = BM["address_rw"];
  SubBMAcc_RW<unsigned char> macc(xpBM, rowInd, colInd, 1);

  size_t n = macc.nrow();
  size_t m = macc.ncol();

  int c, c1, c2;
  std::vector<size_t> ind_na; ind_na.reserve(n);
  unsigned char imputed;

  // tabulate
  for (size_t j = 0; j < m; j++) {

    c = n;
    c1 = c2 = 0;
    ind_na.clear();
    for (size_t i = 0; i < n; i++) {
      unsigned char geno = macc(i, j);
      if (geno == 0) {
        // c0++:
      } else if (geno == 1) {
        c1++;
      } else if (geno == 2) {
        c2++;
      } else {
        ind_na.push_back(i);
        c--;
      }
    }

    // Imputation
    if (ind_na.size() > 0) {

      if (method == 4) { // sample according to allele frequency
        double af = (0.5 * c1 + c2) / c;
        for (auto& k : ind_na) macc(k, j) = ::Rf_rbinom(2, af) + 4;
      } else {
        if (method == 1) { // mode
          imputed = 0;
          int c0 = c - (c1 + c2);
          if (c1 > c0) imputed = 1;
          if (imputed == 0 && c2 > c0) imputed = 2;
          if (imputed == 1 && c2 > c1) imputed = 2;
          imputed += 4;
        } else if (method == 2) { // mean (rounded to 0 decimal place)
          double mean = (c1 + 2.0 * c2) / c;
          imputed = ::Rf_fround(mean, 0) + 4;
        } else if (method == 3) { // mean (rounded to 2 decimal places)
          double mean = (c1 + 2.0 * c2) / c;
          imputed = ::Rf_fround(100 * mean, 0) + 7;
        } else {
          Rcpp::stop("Parameter 'method' should be 1, 2, 3, or 4.");
        }

        for (auto& k : ind_na) macc(k, j) = imputed;
      }
    }

  }

}

/******************************************************************************/
