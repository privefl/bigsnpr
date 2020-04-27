#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List dgTMatrix_to_list(const S4& mat, int offset = 0) {

  IntegerVector ind_row = mat.slot("i");
  IntegerVector ind_col = mat.slot("j");
  NumericVector val     = mat.slot("x");
  int dim = as<IntegerVector>(mat.slot("Dim"))[1];
  List res(dim);

  int K = ind_col.size();
  int cur_ind = 0, last_ind = 0;

  for (int j = 0; j < dim; j++) {

    while (last_ind < K && ind_col[last_ind] == j) last_ind++;

    int n_ind = last_ind - cur_ind;
    IntegerVector i(n_ind);
    NumericVector x(n_ind);

    for (int k = 0; cur_ind < last_ind; k++, cur_ind++) {
      i[k] = ind_row[cur_ind] + offset;
      x[k] = val[cur_ind];
    }

    res[j] = List::create(i, x);
  }

  return res;
}
