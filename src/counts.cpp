// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/isna.hpp>

using namespace Rcpp;


/******************************************************************************/

// [[Rcpp::export]]
ListOf<IntegerMatrix> mycount(SEXP pBigMat,
                              const IntegerVector& indCase,
                              const IntegerVector& indControl) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<char> macc(*xpMat);

  int nCase = indCase.size();
  int nControl = indControl.size();
  int m = xpMat->ncol();

  IntegerMatrix resCase(3, m);
  IntegerMatrix resControl(3, m);

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < nCase; i++) {
      (resCase(macc[j][indCase[i]-1], j))++;
    }
    for (int i = 0; i < nControl; i++) {
      (resControl(macc[j][indControl[i]-1], j))++;
    }
  }

  return(List::create(_["cases"] = resCase,
                      _["controls"] = resControl));
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

  IntegerMatrix resCase(3, m);
  IntegerMatrix resControl(3, m);
  IntegerVector res2(nCase+nControl);

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < nCase; i++) {
      ind = indCase[i]-1;
      tmp = macc[j][ind];
      if (isna(tmp)) {
        (res2[ind])++;
      } else {
        (resCase(tmp, j))++;
      }
    }
    for (int i = 0; i < nControl; i++) {
      ind = indControl[i]-1;
      tmp = macc[j][ind];
      if (isna(tmp)) {
        (res2[ind])++;
      } else {
        (resControl(tmp, j))++;
      }
    }
  }

  return(List::create(_["cols.cases"] = resCase,
                      _["cols.controls"] = resControl,
                      _["rows"] = res2));
}

/******************************************************************************/
