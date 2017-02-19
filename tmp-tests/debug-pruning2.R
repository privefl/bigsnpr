X.cor <- cor(X[, 1:100])
which(X.cor[, 40]^2 > 0.5)
ms <- colMeans(popres$genotypes[, 1:100])
maf <- pmin(ms/2, 1-ms/2)
summary(maf)

i <- 1

ind <- which(X.cor[, i]^2 > 0.5)
ind <- ind[ind > i]
maf[i] < maf[ind]
i <- i + 1

### keep:
# 1, 2, 5-15, 17, 18, 20, 21, 24, 26, 31-35, 38

### remove:
# 1  -> 3, 4
# 17 -> 16
# 17 -> 19
# 20 -> 22
# 21 -> 25, 27, 28, 30
# 24 -> 23
# 24 -> 29
# 34 -> 36, 37
# NA -> 39
