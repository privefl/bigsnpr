set.seed(1)

test <- snp_readExample() # need another example

# generating random phenotypes
print(n <- nrow(test$fam))
print(m <- nrow(test$map))
test$fam$affection <- sample(1:2, size = n, replace = TRUE)

# pruning / clumping with MAF
ind.keep <- snp_pruning(test)
print(length(ind.keep) / m)

stats.MAF <- function(X, ind.train) {
  n <- length(ind.train)
  p <- bigstatsr::big_colstats(X, ind.train)$sum / (2 * n)
  maf <- pmin(p, 1 - p)

  list(S = maf, pS = rank(-maf) / (length(maf) + 1))
} # Here, p-values are used for ranking
X <- test$genotypes
ind <- seq(nrow(X))
maf <- stats.MAF(X, ind)
plot(maf$S, maf$pS, main = "pS is decreasing with respect to S")

print(system.time(
  ind.keep2.1 <- snp_clumping(test, fun.stats = stats.MAF, thr.corr = 0.05)
))
print(system.time(
  ind.keep2.2 <- snp_clumping2(test, S = maf$S, thr.corr = 0.05)
))
all.equal(ind.keep2.1, ind.keep2.2)


print(system.time(
  ind.keep3 <- snp_clumping(test, fun.stats = stats.MAF,
                            thr.pvalue = -log10(0.8))
))
abline(v = min(maf$S[ind.keep3]), h = 0.8, col = "red") # useless?

zcatt <- snp_MAX3(test, val = 0.5)
R2 <- bigstatsr::big_univRegLin(X, test$fam$affection)["R2", ]
plot(zcatt$S, n * R2) # same thing
