#include <Rcpp.h>

using namespace Rcpp;


/******************************************************************************/

// [[Rcpp::export]]
ListOf<IntegerMatrix> mycount(const IntegerMatrix& mat,
                              const IntegerVector& indCase,
                              const IntegerVector& indControl) {
  int nCase = indCase.size();
  int nControl = indControl.size();
  int m = mat.cols();

  // indices begin at 1 in R and 0 in C++
  IntegerVector cases = indCase - 1;
  IntegerVector controls = indControl - 1;

  IntegerMatrix resCase(3, m);
  IntegerMatrix resControl(3, m);

  for (int j = 0; j < m; j++) {
    for (int i = 0; i < nCase; i++) {
      (resCase(mat(cases[i], j), j))++;
    }
    for (int i = 0; i < nControl; i++) {
      (resControl(mat(controls[i], j), j))++;
    }
  }

  return(List::create(_["cases"] = resCase,
                      _["controls"] = resControl));
}
