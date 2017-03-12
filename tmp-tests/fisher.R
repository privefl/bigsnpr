# ASSERT:
# 'x' is 2x2 integer matrix, y = NULL
# hybrid = FALSE, conf.int = FALSE, or = 1, alternative = "two.sided"
fisher <- function (nbNA.ca, nbNA, N.ca, N, relErr = 1 + 1e-7) {

  len <- length(nbNA)
  arr <- array(dim = c(2, 2, len))
  arr[1, 1, ] <- nbNA.ca
  arr[2, 1, ] <- nbNA.co
  arr[1, 2, ] <- N.ca - nbNA.ca
  arr[2, 2, ] <- N.co - nbNA.co

  m <- sum(x[, 1L])
  n <- sum(x[, 2L])
  k <- sum(x[1L, ])
  x <- x[1L, 1L]
  lo <- max(0L, k - n)
  hi <- min(k, m)
  d <- dhyper(lo:hi, m, n, k)
  d.sum <- sum(d)

  relErr <- 1 + 10^(-7)
  PVAL <- sum(d[d <= d[x - lo + 1] * relErr]) / d.sum
}


print(system.time(
  p.fisher2 <- apply(arr, 3, fisher)
)) # 8 sec -> 5 -> 4.5 -> 4
print(all.equal(p.fisher, p.fisher2))



fisher2 <- function(nbNA.ca, nbNA, N.ca, N, relErr = 1 + 1e-7) {
  # precomputation
  dhyper.precomputed <- vector("list", max(nbNA))
  for (m in sort(unique(nbNA))) {
    n <- N - m
    lo <- max(0L, N.ca - n)
    hi <- min(N.ca, m)
    d <- dhyper(lo:hi, m, n, N.ca)
    dhyper.precomputed[[m]] <- d / sum(d)
  }

  # compute p-values
  len <- length(nbNA)
  PVAL <- numeric(len)
  ind <- nbNA.ca - pmax(0L, N.ca - N + nbNA) + 1
  for (i in seq_len(len)) {
    d <- dhyper.precomputed[[nbNA[i]]]
    PVAL[i] <- sum(d[d <= (d[ind[i]] * relErr)])
  }

  PVAL
}

print(system.time(
  p.fisher3 <- fisher2(nbNA.ca, nbNA = nbNA.ca + nbNA.co, N.ca, N)
)) # 8 sec -> 5 -> 4.5 -> 4
print(all.equal(p.fisher, p.fisher3))
