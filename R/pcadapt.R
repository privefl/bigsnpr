################################################################################

getD <- function(X) {
  # robust::covRob(X, estim = "pairwiseGK")$dist
  ogk <- robustbase::covOGK(X, sigmamu = robustbase::s_mad)
  X2 <- sweep(X, 2, ogk$wcenter)
  rowSums((X2 %*% solve(ogk$wcov)) * X2)
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
#' @inherit bigstatsr::linRegPcadapt references
#'
#' @examples
#' test <- snp_attachExtdata()
#' G <- test$genotypes
#' obj.svd <- big_SVD(G, fun.scaling = snp_scaleBinom(), k = 10)
#' plot(obj.svd) # there seems to be 3 "significant" components
#' pcadapt <- snp_pcadapt(G, obj.svd$u[, 1:3])
#' snp_qq(snp_gc(pcadapt))
snp_pcadapt <- function(G, U.row, ind.row = rows_along(G)) {

  K <- ncol(U.row)
  stopifnot(all.equal(crossprod(U.row), diag(K)))

  zscores <- linRegPcadapt(attach.BM(G), U = U.row, rowInd = ind.row)

  fun.pred <- eval(parse(text = sprintf(
    "function(xtr) {
       stats::pchisq(xtr, df = %d, lower.tail = FALSE, log.p = TRUE) / log(10)
     }", K)))

  structure(data.frame(score = getD(zscores)),
            class = c("mhtest", "data.frame"),
            transfo = identity,
            predict = fun.pred)
}

################################################################################
