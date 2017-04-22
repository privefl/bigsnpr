library(bigsnpr)
ex <- snp_attachExtdata()

G <- ex$genotypes
G2 <- attach.BM(G)



corr <- cor(G2[,])
# ind <- which(corr^2 > 0.2, arr.ind = TRUE)
# ind.excl <- ind[ind[, 1] > ind[, 2], ]
# corr[ind.excl]

true <- match(scan("inst/extdata/plink.prune.in", what = ""), ex$map$marker.ID)


ind.row <- rows_along(G)
ind.chr <- cols_along(G)
stats <- big_colstats(G2, ind.row, ind.chr)
m.chr <- length(ind.chr)
keep <- rep(TRUE, m.chr)
n <- length(ind.row)
p <- stats$sum / (2 * n)
maf <- pmin(p, 1 - p)

corr2 <- bigsnpr:::corMat(G2, rowInd = ind.row, colInd = ind.chr,
                          size = 50, thr = rep(0.2, n))

str(obj.svd <- big_SVD(G, big_scale()))
Sys.sleep(2)
str(
  test <- which(bigsnpr:::pruning(G2,
                            rowInd = ind.row,
                            colInd = ind.chr,
                            keep,
                            mafX = maf,
                            sumX = stats$sum,
                            denoX = (n - 1) * stats$var,
                            size = 50,
                            thr = 0.2))
)
Sys.sleep(2)
str(
  test <- bigsnpr:::pruningChr(G,
                               ind.chr = cols_along(G),
                               ind.row = rows_along(G),
                               size = 50,
                               is.size.in.bp = FALSE,
                               infos.pos = NULL,
                               thr.r2 = 0.2,
                               exclude = NULL,
                               nploidy = 2)
)
Sys.sleep(2)
str(ind.keep <- snp_pruning(G, ex$map$chromosome, thr.r2 = 0.2))
Sys.sleep(2)
str(ind.keep <- snp_clumping(G, ex$map$chromosome, thr.r2 = 0.2))



