getD <- function(X) {
  # covRob(X, estim = "pairwiseGK")$dist
  ogk <- robustbase::covOGK(X, sigmamu = robustbase::s_mad)
  X2 <- sweep(X, 2, ogk$wcenter)
  rowSums((X2 %*% solve(ogk$wcov)) * X2)
}

#' Title
#'
#' @param G
#' @param U.row
#' @param ind.row
#'
#' @return
#' @export
#'
#' @examples
snp_pcadapt <- function(G, U.row, ind.row = rows_along(G)) {

  K <- ncol(U.row)
  stopifnot(all.equal(crossprod(U.row), diag(K)))

  zscores <- bigstatsr::linRegPcadapt(attach.BM(G), U = U.row, rowInd = ind.row)

  fun.pred <- eval(parse(text = sprintf(
    "function(xtr) {
       stats::pchisq(xtr, df = %d, lower.tail = FALSE, log.p = TRUE) / log(10)
     }", K)))

  structure(data.frame(score = getD(zscores)),
            class = c("mhtest", "data.frame"),
            transfo = identity,
            predict = fun.pred)
}
