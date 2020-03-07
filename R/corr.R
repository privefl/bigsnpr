################################################################################

#' Correlation
#'
#' Get significant correlations between nearby SNPs of the same chromosome
#' (p-values are computed using a two-sided t-test).
#'
#' @inheritParams bigsnpr-package
#' @param size For one SNP, window size around this SNP to compute correlations.
#' Default is `500`. If not providing `infos.pos` (`NULL`, the default), this is
#' a window in number of SNPs, otherwise it is a window in kb (genetic distance).
#' @param alpha Type-I error for testing correlations.
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
                    fill.diag = TRUE,
                    infos.pos = NULL,
                    ncores = 1) {

  check_args()

  if (is.null(infos.pos)) infos.pos <- 1000 * seq_along(ind.col)
  assert_lengths(infos.pos, ind.col)
  assert_sorted(infos.pos)

  # Get significance thresholds with type-I error `alpha`
  suppressWarnings(
    q.alpha <- stats::qt(p = alpha / 2,
                         df = seq_along(ind.row) - 2,
                         lower.tail = FALSE)
  )
  THR <- q.alpha / sqrt(seq_along(ind.row) - 2 + q.alpha^2)

  corr <- corMat(
    BM     = Gna,
    rowInd = ind.row,
    colInd = ind.col,
    size   = size * 1000,
    thr    = THR,
    pos    = infos.pos,
    ncores = ncores
  )

  corr <- forceSymmetric(corr)
  if (fill.diag) diag(corr) <- 1

  corr
}

################################################################################
