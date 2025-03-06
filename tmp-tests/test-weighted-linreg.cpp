/******************************************************************************/

// [[Rcpp::depends(bigstatsr, rmio)]]

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/BMCodeAcc.h>
#include <bigstatsr/utils.h>

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using std::size_t;

/******************************************************************************/

// [[Rcpp::export]]
ListOf<NumericVector> univLinRegW(Environment BM,
                                  const IntegerVector& rowInd,
                                  const IntegerVector& colInd,
                                  const arma::mat& U,
                                  const arma::vec& y,
                                  const arma::vec& w2,
                                  int ncores = 1) {

  XPtr<FBM> xpBM = BM["address"];
  SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);

  size_t n = macc.nrow();
  size_t m = macc.ncol();
  int K = U.n_cols;
  myassert_size(U.n_rows, n);
  myassert_size(y.n_elem, n);

  arma::vec y2 = w2 % y;
  y2 = y2 - U * (U.t() * y2);
  double y2_sumSq = dot(y2, y2);

  NumericVector betas(m), var(m);

  int chunk_size = ceil(m / (10.0 * ncores));

  #pragma omp parallel num_threads(ncores)
  {
    arma::vec x2(K);

    #pragma omp for schedule(dynamic, chunk_size)
    for (size_t j = 0; j < m; j++) {

      double beta_num = 0, beta_deno = 0;
      x2.zeros();

      for (size_t i = 0; i < n; i++) {
        double x_i = macc(i, j) * w2[i];
        beta_num  += x_i * y2[i];
        beta_deno += x_i * x_i;
        for (int k = 0; k < K; k++) x2[k] += U(i, k) * x_i; // x2 = U.t() * x
      }

      beta_deno -= dot(x2, x2);
      double beta = beta_num / beta_deno;
      double RSS = y2_sumSq - beta_num * beta;

      betas[j] = beta;
      var[j] = RSS / (beta_deno * (n - K - 1));
    }
  }

  return List::create(_["estim"]   = betas,
                      _["std.err"] = sqrt(var));
}

/******************************************************************************/

/*** R
big_univLinReg_w <- function(X, y.train, w.train,
                             ind.train = rows_along(X), ind.col = cols_along(X),
                             covar.train = NULL, thr.eigval = 1e-04, ncores = 1)
{
  bigstatsr:::check_args()
  n <- length(ind.train)
  w2 <- sqrt(w.train)
  covar.train <- sweep(cbind(rep(1, n), covar.train), 1, w2, '*')
  bigassertr::assert_lengths(ind.train, y.train, rows_along(covar.train))
  stopifnot(n > ncol(covar.train))
  SVD <- svd(covar.train, nv = 0)
  eigval.scaled <- SVD$d/(sqrt(n) + sqrt(ncol(covar.train)) - 1)
  K <- sum(eigval.scaled > thr.eigval)
  res <- univLinRegW(BM = X, U = SVD$u[, 1:K, drop = FALSE], w2 = w2,
                     y = y.train, rowInd = ind.train, colInd = ind.col, ncores = ncores)
  res$score <- res$estim/res$std.err
  fun.pred <- eval(parse(text = sprintf("function(xtr) {\n       lpval <- stats::pt(xtr, df = %d, lower.tail = FALSE, log.p = TRUE)\n       (log(2) + lpval) / log(10)\n     }",
                                        n - K - 1)))
  environment(fun.pred) <- baseenv()
  structure(as.data.frame(res), class = c("mhtest", "data.frame"),
            transfo = abs, predict = fun.pred)
}
*/

/******************************************************************************/
