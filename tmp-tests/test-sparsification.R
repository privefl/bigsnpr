library(backbone)
corr0 <- readRDS("tmp-data/corr_chr6_2.rds")
ind <- 1:1000
corr <- abs(corr0)
corr2 <- corr; diag(corr2) <- 0

# corr <- corr2
# Eq (2) of www.pnas.org/cgi/doi/10.1073/pnas.0808904106
alpha <- 0.3
K <- colSums(corr > 0.01)
(thr <- (1 - alpha^(1/(K - 1))) * colSums(corr))
hist(thr, "FD")

cond <- sweep(corr2, 2, thr, '>')
(cond2 <- as((cond + t(cond)) > 0 + 0, "dgCMatrix"))
identical(disparity(corr, alpha = alpha), as(cond2, "symmetricMatrix"))
sum(cond2) / sum(K)

microbenchmark::microbenchmark(
  A = {
    cond <- sweep(corr2, 2, thr, '>')
    as((cond + t(cond)) > 0 + 0, "dgCMatrix")
  },
  B = disparity(corr, alpha = alpha),
  times = 10
)

plot(cond[, 3567], xlim = c(3000, 5000)); points(abs(corr0[, 3567]), pch = 20, col = "red", cex = 0.8)
