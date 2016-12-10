set.seed(1)

test <- snp_readExample() # need another example
X <- test$genotypes
print(n <- nrow(X))
print(m <- ncol(X))

# pruning / clumping with MAF
ind.keep <- snp_pruning(test, thr.corr = 0.1)
print(length(ind.keep) / m)

p <- bigstatsr::big_colstats(X)$sum / (2 * n)
maf <- pmin(p, 1 - p)
ind.keep2 <- snp_clumping(test, S = maf, thr.corr = 0.1)

# clumping for PRS
test$fam$affection <- sample(1:2, size = n, replace = TRUE)
zcatt <- snp_MAX3(test, val = 0.5)
R2 <- bigstatsr::big_univRegLin(X, test$fam$affection)["R2", ]
plot(zcatt$S, n * R2) # same thing
ind.keep3 <- snp_clumping(test, S = zcatt$S, thr.corr = 0.1)
print(length(ind.keep3))
# thresholding though the exclude parameter
ind.keep4 <- snp_clumping(test, S = zcatt$S, thr.corr = 0.1,
                          exclude = which(-log10(zcatt$pS) > 1))
print(length(ind.keep4))
