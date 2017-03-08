require(bigsnpr)

popres <- snp_attach("backingfiles/popres.rds")


X2 <- attach.BM(popres$genotypes)
X2@code <- c(0:2, NA, 0:2, seq(0, 2, by = 0.01), rep(NA, 48))
x <- 1.78899
X2[1, 1] # 2

ALL.RAWS <- as.raw(0:255)
x2 <- round(100 * x) + 8
X2[1, 1] <- ALL.RAWS[x2]
X2[1, 1] # 1.79

X2[1, 1] <- ALL.RAWS[3]
X2[1, 1] # 2

X2@code <- c(0:2, rep(NA_real_, 253))
counts <- big_counts(X2)
tmp <- rownames(counts)[4]

nbNA <- counts[is.na(rownames(counts)), ]
