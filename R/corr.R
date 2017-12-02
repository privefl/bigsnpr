################################################################################

#' Correlation
#'
#' Get significant correlations between nearby SNPs.
#'
#' P-values are computed by a two-sided t-test.
#'
#' @inheritParams bigsnpr-package
#' @param alpha Type-I error for testing correlations.
#' @param size For one SNP, number of SNPs at its left and its right to
#' be tested for being correlated with this particular SNP.
#' @param fill.diag Whether to fill the diagonal with 1s (the default)
#' or to keep it as 0s.
#'
#' @return The correlation matrix. This is a sparse symmetric matrix.
#'
#' @import Matrix
#'
#' @examples
#' test <- snp_attachExtdata()
#'
#' corr <- snp_cor(test$genotypes, ind.col = 1:1000)
#' corr[1:10, 1:10]
#'
#' # Sparsity
#' length(corr@x) / length(corr)
#'
#' @export
snp_cor <- function(Gna,
                    ind.row = rows_along(Gna),
                    ind.col = cols_along(Gna),
                    size = 500,
                    alpha = 0.05,
                    fill.diag = TRUE) {

  check_args()

  # Get significance thresholds with type-I error `alpha`
  suppressWarnings(
    q.alpha <- stats::qt(p = alpha / 2,
                         df = seq_along(ind.row) - 2,
                         lower.tail = FALSE)
  )

  corr <- forceSymmetric(corMat(
    BM = Gna,
    rowInd = ind.row,
    colInd = ind.col,
    size = size,
    thr = q.alpha
  ))

  if (fill.diag) diag(corr) <- 1

  corr
}

################################################################################
