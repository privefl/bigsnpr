require(bigsnpr)

popres <- snp_attach("backingfiles/popresNA.rds")
G <- popres$genotypes
chrs <- popres$map$chromosome
pos <- popres$map$physical.pos
ind.excl <- snp_indLRLDR(chrs, pos)

ind.keep <- snp_clumping(G, chrs, thr.r2 = 0.2, exclude = ind.excl, ncores = 4)
G.svd <- big_randomSVD(G, fun.scaling = snp_scaleBinom(), ind.col = ind.keep,
                       ncores = 4)
zscores <- linRegPcadapt(attach.BM(G), U = G.svd$u, rowInd = rows_along(G))

require(robust)
print(system.time(
  d <- covRob(zscores, estim = "pairwiseGK")$dist
)) # 12 sec

pval <- pchisq(d, df = ncol(zscores), lower.tail = FALSE)
lamGC <- median(d) / qchisq(0.5, df = ncol(zscores))

pvalGC <- pchisq(d / lamGC, df = ncol(zscores), lower.tail = FALSE)

plot(-log10(pvalGC), pch = 19, cex = 0.5, col = chrs)
abline(h = -log10(5e-8), lty = 2)

ind.keep2 <- snp_clumping(G, chrs, thr.r2 = 0.2, exclude = ind.excl,
                         ncores = 4, is.size.in.kb = TRUE, infos.pos = pos)
mean(ind.keep %in% ind.keep2)

ind.keep3 <- snp_pruning(G, chrs, thr.r2 = 0.2, exclude = ind.excl, ncores = 4)
mean(ind.keep3 %in% ind.keep)

ind.keep4 <- snp_pruning(G, chrs, thr.r2 = 0.2, exclude = ind.excl, size = 1000,
                          ncores = 4, is.size.in.kb = TRUE, infos.pos = pos)
mean(ind.keep3 %in% ind.keep4)
