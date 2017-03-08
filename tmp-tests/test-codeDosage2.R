require(bigsnpr)

# snp_readBed("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300.bed",
#             backingfile = "celiac300")

celiac <- snp_attach("backingfiles/celiac300.rds")

counts <- big_counts(celiac$genotypes)
nbNA <- counts[4, ]
nbNA.chr1 <- nbNA[celiac$map$chromosome == 1]

plot(nbNA.chr1, pch = 19, cex = 0.5)
sum(nbNA.chr1 < 5)
hist(nbNA.chr1[nbNA.chr1 < 20])

celiacSub <- subset(celiac, ind.col = which(nbNA.chr1 < 5))
countsSub <- big_counts(celiacSub$genotypes)
firstImpute <- apply(countsSub, 2, which.max) - 1
meanSub <- (countsSub[2, ] + 2 * countsSub[3, ]) / colSums(countsSub[-4, ])
ALL.RAWS <- as.raw(0:255)
firstImpute <- ALL.RAWS[round(meanSub * 100) + 8]

Xsub <- attach.BM(celiacSub$genotypes)
for (i in 1:ncol(Xsub)) {
  indNA <- which(is.na(Xsub[, i]))
  Xsub[indNA, i] <- firstImpute[i]
}
Xsub@code <- c(0:2, NA, 0:2, seq(0, 2, by = 0.01), rep(NA, 48))
plot(Xsub[, 1])

svd <- big_randomSVD(Xsub, fun.scaling = big_scale(), k = 10)
plot(svd$u)
plot(svd$v[, 3], pch = 19, cex = 0.5)
