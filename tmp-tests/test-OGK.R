snp_pcadapt2 <- function(G, U.row, ind.row = rows_along(G)) {

  K <- ncol(U.row)
  stopifnot(all.equal(crossprod(U.row), diag(K)))

  zscores <- linRegPcadapt(attach.BM(G), U = U.row, rowInd = ind.row)

  fun.pred <- eval(parse(text = sprintf(
    "function(xtr) {
       stats::pchisq(xtr, df = %d, lower.tail = FALSE, log.p = TRUE) / log(10)
     }", K)))

  structure(data.frame(score = covRob(zscores, estim = "pairwiseGK")$dist),
            class = c("mhtest", "data.frame"),
            transfo = identity,
            predict = fun.pred)
}

require(bigsnpr)
require(robust)
test <- snp_attachExtdata()
svd <- big_SVD(test$genotypes, snp_scaleBinom(), k = 10)
tmp <- snp_pcadapt(test$genotypes, svd$u[, 1:3])
tmp2 <- snp_pcadapt2(test$genotypes, svd$u[, 1:3])

tmp3 <- pcadapt::pcadapt(t(attach.BM(test$genotypes)[,]), K = 3)
tmp4 <- tmp3$chi2.stat * tmp3$gif
