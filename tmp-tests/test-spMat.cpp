#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::export]]
sp_mat test_sp() {

  // batch insertion of two values at (5, 6) and (9, 9)
  umat locations;
  locations << 5 << 9 << endr
            << 6 << 9 << endr;

  vec values;
  values << 1.5 << 3.2 << endr;

  sp_mat X(locations, values);


  // batch insertion of two values at (5, 6) and (9, 9)
  umat locations2;
  locations2 << 5 << 9 << endr
             << 7 << 3 << endr;

  vec values2;
  values2 << 10.5 << 30.2 << endr;
  X.elem(locations2) = values2; // DOES NOT WORK

  return X;
}


/*** R
test_sp()
*/
