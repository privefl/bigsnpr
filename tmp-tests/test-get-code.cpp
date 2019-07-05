#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix get_code(int NA_VAL = 3) {

  IntegerVector num = IntegerVector::create(2, NA_VAL, 1, 0);
  IntegerMatrix code(4, 256);

  int i, k, k2;
  int coeff = 1;
  for (i = 0; i < 4; i++) {
    for (k = 0; k < 256; k++) {
      k2 = k / coeff;
      code(i, k) = num[k2 % 4];
    }
    coeff *= 4;
  }

  return code;
}


/*** R
test <- get_code()
true <- bigsnpr:::getCode()
storage.mode(true) <- "integer"
identical(test, true)
*/
