source("R/utils.R")

require(bigsnpr)

celiac <- AttachBigSNP("celiac")
X <- celiac$genotypes
n <- nrow(X)

celiacRaw <- AttachBigSNP("celiacRaw")

x <- list(X = celiacRaw$genotypes, tab = getCode(),
          q = as.integer(seq(0, n-1) %/% 4),
          r = as.integer(seq(0, n-1) %% 4))
class(x) <- "bigRaw"

`[.bigRaw` <- function(x, ind_i = seq(nrow(X)), ind_j) {
  rawToBigPart2((x$X)@address, ind_i, ind_j, x$tab, x$q, x$r)
}


require(microbenchmark)
print(microbenchmark(
  x[1:10, 1:5 + 10000],  # 13 ms
  X[1:10, 1:5 + 10000]   # 54 ms
))

