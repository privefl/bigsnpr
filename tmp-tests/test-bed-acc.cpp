#include "../src/bed-acc.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void test_acc(const std::string path,
              int n_total, int m_total,
              const IntegerVector& row_ind,
              const IntegerVector& col_ind,
              const RawMatrix& lookup_byte) {

  bedAcc macc(path, n_total, m_total, row_ind, col_ind, lookup_byte);
  size_t n = macc.nrow();
  size_t m = macc.ncol();

  Rcout << n << " / " << m << std::endl;
  Rcout << n << " / " << m << " / " << int(macc(0, 0)) << std::endl;
}

/*** R
bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
bimfile <- sub("\\.bed$", ".bim", bedfile)
famfile <- sub("\\.bed$", ".fam", bedfile)
n_total <- bigreadr::nlines(famfile)
m_total <- bigreadr::nlines(bimfile)
ind.row <- seq_len(n_total)
ind.col <- seq_len(m_total)
lookup_byte <- bigsnpr:::getCode()
test_acc(bedfile, n_total, m_total, ind.row, ind.col, lookup_byte)
*/
