################################################################################

context("SCT")

################################################################################

snp <- snp_attachExtdata()
G <- snp$genotypes
CHR <- rep(1:2, c(2542, 2000))
POS <- snp$map$physical.pos

lpval <- -log10(runif(ncol(G)))
betas <- rnorm(ncol(G), sd = 0.1)
y <- sample(0:1, size = nrow(G), replace = TRUE)

################################################################################

all_keep <- snp_grid_clumping(G, CHR, POS, lpS = lpval, ncores = 2,
                              grid.thr.r2 = c(0.05, 0.2, 0.8),
                              grid.base.size = c(100, 200))
expect_length(all_keep, 2)
expect_length(sample(all_keep, 1)[[1]], 6)

grid <- attr(all_keep, "grid")[1:2]
all_keep2 <- lapply(rows_along(grid), function(i) {
  expect_identical(snp_clumping(G, CHR, S = lpval, thr.r2 = grid$thr.r2[i],
                                size = grid$size[i], infos.pos = POS, ncores = 2),
                   c(all_keep[[1]][[i]], all_keep[[2]][[i]]))
})

################################################################################

multi_PRS <- snp_grid_PRS(G, all_keep, betas, lpval, grid.lpS.thr = 0:5, ncores = 2)

multi_PRS2 <- lapply(all_keep2, function(ind.keep) {
  snp_PRS(G, betas[ind.keep], ind.keep = ind.keep, lpS.keep = lpval[ind.keep],
          thr.list = 0:5)
})
expect_equal(do.call("cbind", multi_PRS2), multi_PRS, check.attributes = FALSE)

multi_PRS3 <- lapply(unlist(all_keep, recursive = FALSE), function(ind.keep) {
  snp_PRS(G, betas[ind.keep], ind.keep = ind.keep, lpS.keep = lpval[ind.keep],
          thr.list = 0:5)
})
expect_equal(do.call("cbind", multi_PRS3), big_attach(attr(multi_PRS, "rds"))[],
             check.attributes = FALSE)

################################################################################

new_betas <- snp_grid_stacking(multi_PRS, y, alphas = 1e-3, ncores = 2)
expect_length(new_betas$beta.covar, 0)
expect_equal(
  predict(new_betas$mod, big_attach(attr(multi_PRS, "rds")), proba = FALSE),
  new_betas$intercept + big_prodVec(G, new_betas$beta.G),
  check.attributes = FALSE)

################################################################################
