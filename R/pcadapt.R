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

plot.big_SVD <- function(x, type = c("screeplot", "scores", "loadings"),
                         nval = length(x$d),
                         scores = c(1, 2),
                         loading = 1,
                         pch = 19, cex = 0.5,
                         ...) {

  type <- match.arg(type)
  if (type == "screeplot") {
    plot(x$d[seq_len(nval)], type = "b", ylab = "singular value",
         pch = pch, cex = cex, ...)
  } else if (type == "scores") {
    plot(predict(x)[, scores],
         xlab = paste0("PC", scores[1]),
         ylab = paste0("PC", scores[2]),
         pch = pch, cex = cex, ...)
  } else if (type == "loadings") {
    plot(x$v[, loading], ylab = paste0("loadings of PC", loading),
         pch = pch, cex = cex, ...)
  }
}
plot(svd, nval = 10)
plot(svd, type = "scores", scores = 3:4)
plot(svd, type = "loadings", loading = 1)

