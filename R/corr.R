################################################################################

cor0 <- function(Gna,
                 ind.row = rows_along(Gna),
                 ind.col = cols_along(Gna),
                 size = 500,
                 alpha = 1,
                 thr_r2 = 0,
                 fill.diag = TRUE,
                 infos.pos = NULL,
                 ncores = 1) {

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

  ind_val <- corMat(
    obj    = Gna,
    rowInd = ind.row,
    colInd = ind.col,
    size   = size * 1000,
    thr    = pmax(THR, sqrt(thr_r2)),
    pos    = infos.pos,
    ncores = ncores,
    fill_diag = fill.diag
  )

  m <- length(ind.col)
  # corr <- Matrix::sparseMatrix(
  #   i = unlist(lapply(ind_val, function(.) .$i)),
  #   j = rep(1:m, times = sapply(ind_val, function(.) length(.$i))),
  #   x = unlist(lapply(ind_val, function(.) .$x)),
  #   dims = c(m, m),
  #   symmetric = TRUE)
  corr <- new("dsCMatrix", uplo = "U")
  corr@Dim <- c(m, m)
  corr@i <- unlist(lapply(ind_val, function(.) .$i))
  corr@p <- c(0L, cumsum(sapply(ind_val, function(.) length(.$i))))
  corr@x <- unlist(lapply(ind_val, function(.) .$x))

  # if (fill.diag) diag(corr) <- 1
  # if (corr@uplo == "L" && Matrix::isDiagonal(corr))
  #   corr@uplo <- "U"

  if (anyNA(corr@x))
    warning2("NA or NaN values in the resulting correlation matrix.")

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
#'   Default is `500`. If not providing `infos.pos` (`NULL`, the default), this is
#'   a window in number of SNPs, otherwise it is a window in kb (physical distance).
#'   In case you provide `infos.pos` in centimorgans (genetic distance),
#'   you should divide this by 1000 because it is internally multiplied by 1000
#'   (i.e. use `3 / 1000` for 3 cM).
#' @param alpha Type-I error for testing correlations.
#'   Default is `1` (no threshold is applied).
#' @param thr_r2 Threshold to apply on squared correlations. Default is `0`.
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
                    ncores = 1) {

  args <- as.list(environment())
  names(args)[names(args) == "obj.bed"] <- "Gna"

  check_args()

  do.call(cor0, args)
}

################################################################################
