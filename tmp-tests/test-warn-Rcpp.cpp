#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void test_warn() {
  // Rcpp::warning("%d This is wrong %.", 2);
  Rcpp::warning("%d This is wrong %%.", 2);
}


/*** R
test_warn()
*/
