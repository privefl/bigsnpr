// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/isna.hpp>

using namespace Rcpp;


/******************************************************************************/

// [[Rcpp::export]]
IntegerMatrix mycount(const SEXP pBigMat,
                      const IntegerVector& indCase,
                      const IntegerVector& indControl) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int nCase = indCase.size();
  int nControl = indControl.size();
  int m = xpMat->ncol();

  IntegerMatrix res(6, m);

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < nCase; i++) {
      (res(macc[j][indCase[i]-1], j))++;
    }
    for (int i = 0; i < nControl; i++) {
      (res(macc[j][indControl[i]-1]+3, j))++;
    }
  }

  return(res);
}

/******************************************************************************/

// [[Rcpp::export]]
ListOf<SEXP> mycount2(SEXP pBigMat,
                      const IntegerVector& indCase,
                      const IntegerVector& indControl) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int nCase = indCase.size();
  int nControl = indControl.size();
  int m = xpMat->ncol();

  char tmp;
  int ind;

  IntegerMatrix res(6, m);
  IntegerVector res2(nCase+nControl);

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < nCase; i++) {
      ind = indCase[i]-1;
      tmp = macc[j][ind];
      if (isna(tmp)) {
        (res2[ind])++;
      } else {
        (res(tmp, j))++;
      }
    }
    for (int i = 0; i < nControl; i++) {
      ind = indControl[i]-1;
      tmp = macc[j][ind];
      if (isna(tmp)) {
        (res2[ind])++;
      } else {
        (res(tmp+3, j))++;
      }
    }
  }

  return(List::create(_["counts.col"] = res,
                      _["counts.row"] = res2));
}

/******************************************************************************/
