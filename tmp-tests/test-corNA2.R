rowMeans(replicate(200, {
  X <- sapply(runif(20, min = 0.05, max = 0.5), function(af) {
    rbinom(1000, size = 2, prob = af)
  })

  corr0 <- cor(X)

  indNA <- sample(length(X), length(X) * 0.2)
  XNA <- X; XNA[indNA] <- NA

  corr <- cor(XNA, use = "pairwise.complete.obs")
  XNA_scaled <- scale(XNA); XNA_scaled <- ifelse(is.na(XNA_scaled), 0, XNA_scaled)
  corr2 <- cor(XNA_scaled)

  corr[1:5, 1:5]
  corr2[1:5, 1:5]
  corr3 <- crossprod(XNA_scaled)
  corr4 <- corr3 / tcrossprod(sqrt(diag(corr3)))
  corr4[1:5, 1:5]
  stopifnot(isTRUE(all.equal(corr2, corr4)))
  corr5 <- crossprod(sweep(XNA_scaled, 2, sqrt(colSums(!is.na(XNA)) - 1), '/'))
  corr5[1:5, 1:5]
  stopifnot(isTRUE(all.equal(corr5, corr4)))

  c(
    readr::parse_number(all.equal(corr,  corr0)),
    # readr::parse_number(all.equal(corr3, corr0)),
    readr::parse_number(all.equal(corr4, corr0))
  )
}))
