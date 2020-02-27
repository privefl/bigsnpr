################################################################################

part_mult_lin_reg <- function(X, ind, ind.row, U) {
  multLinReg(X, ind_row = ind.row, ind_col = ind, U = U)
}

pcadapt0 <- function(G, U.row, ind.row, ind.col, ncores) {

  if (is.null(dim(U.row))) U.row <- as.matrix(U.row)  # vector
  K <- ncol(U.row)
  assert_lengths(rows_along(U.row), ind.row)

  tscores <- big_parallelize(G, p.FUN = part_mult_lin_reg,
                             p.combine = "rbind", ncores = ncores,
                             ind = ind.col, ind.row = ind.row, U = U.row)

  dist <- `if`(K == 1, (drop(tscores) - stats::median(tscores))^2,
               bigutilsr::dist_ogk(tscores))

  fun.pred <- eval(parse(text = sprintf(
    "function(xtr) {
       stats::pchisq(xtr, df = %d, lower.tail = FALSE, log.p = TRUE) / log(10)
     }", K)))
  environment(fun.pred) <- baseenv()

  snp_gc(
    structure(data.frame(score = dist),
              class = c("mhtest", "data.frame"),
              transfo = identity,
              predict = fun.pred)
  )
}

################################################################################

#' Outlier detection
#'
#' Method to detect genetic markers involved in biological adaptation.
#' This provides a statistical tool for outlier detection based on
#' Principal Component Analysis. This corresponds to the statistic based
#' on mahalanobis distance, as implemented in package **pcadapt**.
#'
#' @inheritParams bigsnpr-package
#' @param U.row Left singular vectors (not scores, \eqn{U^T U = I})
#' corresponding to `ind.row`.
#'
#' @return An object of classes `mhtest` and `data.frame` returning one
#' score by SNP. See `methods(class = "mhtest")`.
#' @seealso [snp_manhattan], [snp_qq] and [snp_gc].
#' @export
#'
#' @references
#' Luu, K., Bazin, E., & Blum, M. G. (2017).
#' pcadapt: an R package to perform genome scans for selection
#' based on principal component analysis.
#' Molecular ecology resources, 17(1), 67-77.
#'
#' @examples
#' test <- snp_attachExtdata()
#' G <- test$genotypes
#' obj.svd <- big_SVD(G, fun.scaling = snp_scaleBinom(), k = 10)
#' plot(obj.svd) # there seems to be 3 "significant" components
#' pcadapt <- snp_pcadapt(G, obj.svd$u[, 1:3])
#' snp_qq(snp_gc(pcadapt))
#'
snp_pcadapt <- function(G, U.row,
                        ind.row = rows_along(G),
                        ind.col = cols_along(G),
                        ncores = 1) {
  check_args()
  pcadapt0(G, U.row, ind.row, ind.col, ncores)
}

################################################################################

#' @export
#' @rdname snp_pcadapt
bed_pcadapt <- function(obj.bed, U.row,
                        ind.row = rows_along(obj.bed),
                        ind.col = cols_along(obj.bed),
                        ncores = 1) {
  check_args()
  pcadapt0(obj.bed$light, U.row, ind.row, ind.col, ncores)
}

################################################################################
