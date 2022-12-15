#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List test_xgboost(const List& params, const NumericMatrix& X, const NumericVector& y) {

  // Obtaining namespace of package
  // Environment pkg = ;

  // Picking up function from package
  Function xgb_cv = Environment::namespace_env("xgboost")["xgb.cv"];

  // Executing function
  return xgb_cv(params, X, 20, 10, y, NA_REAL, true);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
library(xgboost)

data(agaricus.train, package='xgboost')
data(agaricus.test, package='xgboost')
train <- agaricus.train
test <- agaricus.test

mod <- xgboost::xgb.cv(list(), train$data, 20, 10, train$label, prediction = TRUE)
mod$pred
plot(mod$pred, train$label)
mod2 <- test_xgboost(list(), as.matrix(train$data), train$label)
predict(mod2, as.matrix(test$data))
*/
