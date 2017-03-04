# sqrt(n) + sqrt(m) - 1

require(bigsnpr)

popres <- snp_attach("backingfiles/popres.desc")

ind <- snp_indLRLDR(popres)
maf <- snp_MAF(popres$genotypes)
ind.keep <- snp_clumping(popres, S = maf, thr.r2 = 0.2,
                         exclude = union(ind, which(maf < 0.05)), ncores = 3)

X.svd <- big_randomSVD(popres$genotypes, fun.scaling = big_scale(),
                       ind.col = ind.keep, k = 200, verbose = TRUE, ncores = 3)


THR <- sqrt(nrow(X.svd$u)) + sqrt(nrow(X.svd$v)) - 1
plot(X.svd$d); abline(h = THR, col = "red")

ind2 <- which(X.svd$d > THR)

X2 <- deepcopy(popres$genotypes, cols = ind.keep)
test2 <- svd(scale(X2[,]), nu = 200, nv = 200)
