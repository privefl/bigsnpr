library(bigsnpr)
snp <- snp_attachExtdata()
G <- snp$genotypes
G2 <- big_copy(G, ind.row = rep(rows_along(G), 10),
               ind.col = rep(cols_along(G), 1))
dim(G2)

set.seed(1); ord <- sample(ncol(G2))
pos <- cols_along(G2)
stats <- big_colstats(G2)

RcppParallel::setThreadOptions(1)
system.time(
  ind.keep0 <- bigsnpr:::clumping_chr(
    G2, rows_along(G2), cols_along(G2), ord, pos,
    stats$sum, (nrow(G2) - 1) * stats$var, 500, 0.2)
) # 24 sec

for (ncores in c(1, 2, 4)) {
  RcppParallel::setThreadOptions(ncores)
  Rcpp::sourceCpp('tmp-tests/test-clumping-parallel.cpp')
  remain <- rep(TRUE, ncol(G2))
  print(system.time(
    ind.keep <- clumping(G2, rows_along(G2), cols_along(G2), ord, pos, remain,
                         stats$sum, (nrow(G2) - 1) * stats$var, 500, 0.2)
  ))
  stopifnot(identical(ind.keep, ind.keep0))
}

## 0.01:
# 1 -> 1.5 sec
# 2 -> 0.7 sec
# 4 -> 0.4 sec
## 0.2:
# 1 -> 23 sec
# 2 -> 12 sec
# 4 ->  5 sec
## 0.9:
# 1 -> 26
# 2 -> 16
# 4 -> 6
