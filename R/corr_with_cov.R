################################################################################

cor_with_cov0 <- function(G, Z,
                          ind.row = rows_along(G),
                          ind.col = cols_along(G),
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

  ind_val <- corMatCov(
    obj    = G,
    rowInd = ind.row,
    colInd = ind.col,
    Z      = Z,
    size   = size * 1000,
    thr    = pmax(THR, sqrt(thr_r2)),
    pos    = infos.pos,
    ncores = ncores,
    fill_diag = fill.diag
  )

  m <- length(ind.col)
  corr <- new("dsCMatrix", uplo = "U")
  corr@Dim <- c(m, m)
  corr@i <- unlist(lapply(ind_val, function(.) .$i))
  corr@p <- c(0L, cumsum(sapply(ind_val, function(.) length(.$i))))
  corr@x <- unlist(lapply(ind_val, function(.) .$x))

  if (anyNA(corr@x))
    warning2("NA or NaN values in the resulting correlation matrix.")

  corr
}

################################################################################

#' @param covar.row Matrix of covariables (corresponding to `ind.row`) to adjust
#'   for when computing correlations. You can use [bigstatsr::covar_from_df()]
#'   to convert from a data frame.
#'
#' @rdname snp_cor
#' @export
#'
snp_cor_with_cov <- function(G, covar.row,
                             ind.row = rows_along(G),
                             ind.col = cols_along(G),
                             size = 500,
                             alpha = 1,
                             thr_r2 = 0,
                             fill.diag = TRUE,
                             infos.pos = NULL,
                             ncores = 1) {

  args <- as.list(environment())

  check_args()
  assert_lengths(rows_along(covar.row), ind.row)

  # get SVD of covar (+ intercept)
  SVD <- svd(cbind(1, covar.row), nv = 0)
  eigval.scaled <- SVD$d / (sqrt(nrow(covar.row)) + sqrt(ncol(covar.row)) - 1)
  U <- SVD$u[, eigval.scaled > 1e-4, drop = FALSE]

  args[["Z"]] <- t(big_cprodMat(G, U, ind.row, ind.col))
  args[["covar.row"]] <- NULL

  do.call(cor_with_cov0, args)
}

################################################################################

#' #' @rdname snp_cor
#' #' @export
#' bed_cor <- function(obj.bed,
#'                     ind.row = rows_along(obj.bed),
#'                     ind.col = cols_along(obj.bed),
#'                     size = 500,
#'                     alpha = 1,
#'                     thr_r2 = 0,
#'                     fill.diag = TRUE,
#'                     infos.pos = NULL,
#'                     ncores = 1) {
#'
#'   args <- as.list(environment())
#'   names(args)[names(args) == "obj.bed"] <- "G"
#'
#'   check_args()
#'
#'   do.call(cor0, args)
#' }

################################################################################
