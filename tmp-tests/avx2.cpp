#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
bool has_avx2() {
  // https://stackoverflow.com/a/25833114/6103040
  return __builtin_cpu_supports("avx2");
}

/*** R
has_avx2()
*/
