// [[Rcpp::export]]
void computeK(mat& source, NumericVector norms, double lambda) {
  int nrows = source.n_rows;
  int ncols = source.n_cols;

  for (int j = 0; j < ncols; j++) {
    for (int i = 0; i < nrows; i++) {
      source(i,j) *= 2;
      source(i,j) -= norms[i] + norms[j];
    }
  }

  source = exp(lambda * source);

  return;
}

// [[Rcpp::export]]
void computeK2(mat& source, NumericVector norms) {
  int nrows = source.n_rows;
  int ncols = source.n_cols;

  for (int j = 0; j < ncols; j++) {
    for (int i = 0; i < nrows; i++) {
      source(i,j) *= 2;
      source(i,j) -= norms[i] + norms[j];
    }
  }

  return;
}

// [[Rcpp::export]]
NumericVector& mloss2(const mat& source, const NumericVector& y,
                      const NumericVector& lambdas, NumericVector& sums) {

  int nrows = source.n_rows;
  int ncols = source.n_cols;
  int nlambdas = lambdas.size();
  double tmp;

  for (int j = 0; j < ncols; j++) {
    for (int i = j+1; i < nrows; i++) {
      if (y[i] != y[j]) {
        for (int k = 0; k < nlambdas; k++) {
          sums[k] += exp(2 * lambdas[k] * source(i,j));
        }
      } else {
        for (int k = 0; k < nlambdas; k++) {
          tmp = 1 - exp(lambdas[k] * source(i,j));
          sums[k] += tmp * tmp;
        }
      }
    }
  }

  return(sums);
}


// [[Rcpp::export]]
mat& computeK3(mat& source, const NumericVector& normsRow,
               const NumericVector& normsCol, double lambda) {
  int nrows = source.n_rows;
  int ncols = source.n_cols;

  for (int j = 0; j < ncols; j++) {
    for (int i = 0; i < nrows; i++) {
      source(i,j) *= 2;
      source(i,j) -= normsRow[i] + normsCol[j];
    }
  }

  source = exp(lambda * source);

  return(source);
}



// [[Rcpp::export]]
void toKernel(SEXP pBigMat, const NumericVector& norms, double lambda) {
  XPtr<BigMatrix> xpMat(pBigMat);
  MatrixAccessor<double> macc(*xpMat);

  int n = xpMat->nrow();

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      macc[j][i] = exp(lambda * (2.0 * macc[j][i] - norms[i] - norms[j]));
    }
  }

  return;
}

