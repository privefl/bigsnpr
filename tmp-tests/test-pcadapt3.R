snp_pcadapt <- function(G, U.row, ind.row = rows_along(G)) {

  K <- ncol(U.row)
  stopifnot(all.equal(crossprod(U.row), diag(K)))

  zscores <- linRegPcadapt(attach.BM(G), U = U.row, rowInd = ind.row)
  d <- covRob(zscores, estim = "pairwiseGK")$dist

  fun.pred <- eval(parse(text = sprintf(
    "function(xtr) stats::pchisq(xtr, df = %d, lower.tail = FALSE)", K)))

  structure(data.frame(score = d),
            class = c("mhtest", "data.frame"),
            transfo = identity,
            predict = fun.pred)
}

tmp <- snp_pcadapt(G, G.svd$u)
snp_qq(tmp)
snp_qq(snp_gc(tmp))
snp_manhattan(snp_gc(tmp), popres$map)

plot(G.svd$d, type = "b")

tmp <- snp_pcadapt(G, G.svd$u[, 1:5])
snp_qq(tmp)
snp_qq(snp_gc(tmp))
snp_manhattan(snp_gc(tmp), popres$map)
