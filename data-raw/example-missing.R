library(bigsnpr)

gen <- function(n, m) {
  I <- 1:m
  p <- I / (2 * m + 1)

  mat <- outer(I, I, FUN = function(i, j) {
    1 / (abs(i - j) + 1)
  })

  bindata::rmvbin(n, p, bincorr = mat) +
    bindata::rmvbin(n, p, bincorr = mat)
}

N <- 200
M <- 500

fake <- snp_fake(N, M)
G <- fake$genotypes
G[] <- rep(gen(N, M / 20), 20)

nbNA <- VGAM::rbetabinom.ab(M, size = N, shape1 = 0.6, shape2 = 20)
indNA <- cbind(
  unlist(
    lapply(nbNA, function(nb) {
      `if`(nb > 0, sample(N, size = nb), NULL)
    })
  ),
  rep(cols_along(G), nbNA)
)
G[indNA] <- as.raw(3)

fake$map$chromosome <- sort(sample(1:2, M, TRUE))
snp_writeBed(fake, "inst/extdata/example-missing.bed")
