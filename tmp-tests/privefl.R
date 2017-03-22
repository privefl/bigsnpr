require(bigsnpr)
require(pcadapt)

test <- snp_attachExtdata()
mat <- t(attach.BM(test$genotypes)[,])

svd <- big_SVD(test$genotypes, snp_scaleBinom(), k = 10)

x <- pcadapt(mat, K = 10)
plot(x,option="scores")
plot(svd)


x <- pcadapt(mat, K = 3)
tmp <- snp_pcadapt(test$genotypes, svd$u[, 1:3])
plot(x$stat, tmp[[1]])
plot(x$pvalues, predict(snp_gc(tmp)), log = "xy")
abline(0, 1, col = "red")

snp_qq(snp_gc(tmp))
plot(x, option = "qqplot")
plot(x, option = "manhattan")
snp_manhattan(snp_gc(tmp), test$map)
