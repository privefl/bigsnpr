################################################################################

cor0 <- function(Gna,
                 ind.row = rows_along(Gna),
                 ind.col = cols_along(Gna),
                 size = 500,
                 alpha = 1,
                 thr_r2 = 0,
                 fill.diag = TRUE,
                 infos.pos = NULL,
                 info = rep(1, length(ind.col)),
                 ncores = 1) {

  if (is.null(infos.pos)) infos.pos <- 1000 * seq_along(ind.col)
  assert_lengths(infos.pos, info, ind.col)
  assert_sorted(infos.pos)
  if (!all(0.1 <= info & info <= 1))
    stop2("All values of 'info' must be between 0.1 and 1.")

  # Get significance thresholds with type-I error `alpha`
  suppressWarnings(
    q.alpha <- stats::qt(p = alpha / 2,
                         df = seq_along(ind.row) - 2,
                         lower.tail = FALSE)
  )
  THR <- q.alpha / sqrt(seq_along(ind.row) - 2 + q.alpha^2)

  ind_val <- corMat(
    obj    = Gna,
    rowInd = ind.row,
    colInd = ind.col,
    size   = size * 1000,
    thr    = pmax(THR, sqrt(thr_r2)),
    pos    = infos.pos,
    info   = info,
    ncores = ncores
  )

  m <- length(ind.col)
  corr <- Matrix::sparseMatrix(
    i = unlist(lapply(ind_val, function(.) .$i)),
    j = rep(1:m, times = sapply(ind_val, function(.) length(.$i))),
    x = unlist(lapply(ind_val, function(.) .$x)),
    dims = c(m, m),
    symmetric = TRUE)

  if (fill.diag) diag(corr) <- 1
  if (corr@uplo == "L" && Matrix::isDiagonal(corr))
    corr@uplo <- "U"

  corr
}

################################################################################

#' Correlation matrix
#'
#' Get significant (Pearson) correlations between nearby SNPs of the same chromosome
#' (p-values are computed using a two-sided t-test).
#'
#' @inheritParams bigsnpr-package
#' @param size For one SNP, window size around this SNP to compute correlations.
#' Default is `500`. If not providing `infos.pos` (`NULL`, the default), this is
#' a window in number of SNPs, otherwise it is a window in kb (genetic distance).
#' @param alpha Type-I error for testing correlations.
#'   Default is `1` (no threshold is applied).
#' @param thr_r2 Threshold to apply on squared correlations. Default is `0`.
#' @param info Vector of imputation INFO scores to correct correlations when
#'   they are computed from imputed dosage data.
#' @param fill.diag Whether to fill the diagonal with 1s (the default)
#' or to keep it as 0s.
#'
#' @return The (Pearson) correlation matrix. This is a sparse symmetric matrix.
#'
#' @import Matrix
#'
#' @examples
#' test <- snp_attachExtdata()
#' G <- test$genotypes
#'
#' corr <- snp_cor(G, ind.col = 1:1000)
#' corr[1:10, 1:10]
#'
#' # Sparsity
#' length(corr@x) / length(corr)
#'
#' @export
#'
snp_cor <- function(Gna,
                    ind.row = rows_along(Gna),
                    ind.col = cols_along(Gna),
                    size = 500,
                    alpha = 1,
                    thr_r2 = 0,
                    fill.diag = TRUE,
                    infos.pos = NULL,
                    info = rep(1, length(ind.col)),
                    ncores = 1) {

  args <- as.list(environment())

  check_args()

  do.call(cor0, args)
}

################################################################################

#' @rdname snp_cor
#' @export
bed_cor <- function(obj.bed,
                    ind.row = rows_along(obj.bed),
                    ind.col = cols_along(obj.bed),
                    size = 500,
                    alpha = 1,
                    thr_r2 = 0,
                    fill.diag = TRUE,
                    infos.pos = NULL,
                    info = rep(1, length(ind.col)),
                    ncores = 1) {

  args <- as.list(environment())
  names(args)[names(args) == "obj.bed"] <- "Gna"

  check_args()

  do.call(cor0, args)
}

################################################################################
