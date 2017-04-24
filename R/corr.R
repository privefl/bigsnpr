################################################################################

#' Correlation
#'
#' Get significant correlations between nearby SNPs.
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
#' corr <- snp_cor(test$genotypes)
#' corr[1:10, 1:10]
#'
#' length(corr@x) / length(corr)
#'
#' @export
snp_cor <- function(G,
                    ind.row = rows_along(G),
                    ind.col = cols_along(G),
                    size = 500,
                    alpha = 0.05,
                    fill.diag = TRUE) {

  X <- attach.BM(G)

  # get the significance threshold with type-I error `alpha`
  q.alpha <- stats::qchisq(alpha, df = 1, lower.tail = FALSE)

  corr <- forceSymmetric(corMat(
    BM = X,
    rowInd = ind.row,
    colInd = ind.col,
    size = size,
    thr = q.alpha / seq_along(ind.row)
  ))

  if (fill.diag) diag(corr) <- 1

  corr
}

################################################################################
