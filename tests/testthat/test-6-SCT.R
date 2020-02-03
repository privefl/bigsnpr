################################################################################

context("SCT")

################################################################################

test_that("seq_log() works", {
  expect_equal(seq_log(1, 1000, 4), 10^(0:3))
  expect_equal(seq_log(1, 100,  5), 10^(0:4 / 2))
  expect_equal(seq_log(1000, 1, 4), rev(seq_log(1, 1000, 4)))
  expect_equal(seq_log(100,  1, 5), rev(seq_log(1, 100,  5)))
  expect_equal(seq_log(1, 1, 1), 1)
  expect_equal(seq_log(1, 1, 5), rep(1, 5))
  expect_error(seq_log(1, 1000, -4), "'length.out' must be a non-negative number")
})

################################################################################

skip_if(is_cran)

snp <- snp_attachExtdata()
G <- snp$genotypes
CHR <- rep(1:2, c(2542, 2000))
POS <- snp$map$physical.pos

lpval <- -log10(runif(ncol(G)))
betas <- rnorm(ncol(G), sd = 0.1)
y <- sample(0:1, size = nrow(G), replace = TRUE)

################################################################################

expect_error(snp_grid_clumping(G, CHR, sample(POS), lpS = lpval),
             "'pos.chr' is not sorted.")
expect_length(snp_grid_clumping(G, c(CHR[-1], 22), c(POS[-1], 1), lpS = lpval,
                                grid.thr.r2 = 0.2, grid.base.size = 50), 3)

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

infos <- runif(ncol(G), 0.2)
all_keep3 <- snp_grid_clumping(G, CHR, POS, lpS = lpval, ncores = 2,
                               grid.thr.r2 = c(0.05, 0.2, 0.8),
                               grid.base.size = c(100, 200),
                               infos.imp = infos,
                               grid.thr.imp = c(0.3, 0.8, 0.95))
grid3 <- attr(all_keep3, "grid")
expect_equal(dim(grid3), c(18, 4))
expect_equal(grid3$thr.imp, rep(c(0.3, 0.8, 0.95), each = 6))
expect_equal(grid3$grp.num, rep(1, 18))

groups <- lapply(c(0.3, 0.8, 0.95), function(thr) which(infos >= thr))
all_keep4 <- snp_grid_clumping(G, CHR, POS, lpS = lpval, ncores = 2,
                               grid.thr.r2 = c(0.05, 0.2, 0.8),
                               grid.base.size = c(100, 200),
                               groups = groups)
expect_equal(all_keep4, all_keep3, check.attributes = FALSE)
grid4 <- attr(all_keep4, "grid")
expect_equal(dim(grid4), c(18, 4))
expect_equal(grid4$thr.imp, rep(1, 18))
expect_equal(grid4$grp.num, rep(1:3, each = 6))

groups2 <- list(NULL, 1, cols_along(G))
all_keep5 <- snp_grid_clumping(G, CHR, POS, lpS = lpval, ncores = 2,
                               grid.thr.r2 = c(0.05, 0.2, 0.8),
                               grid.base.size = c(100, 200),
                               groups = groups2)
expect_equal(all_keep5, check.attributes = FALSE,
             list(c(rep(list(integer(), 1), each = 6), all_keep[[1]]),
                  c(rep(list(integer()), 12), all_keep[[2]])))
grid5 <- attr(all_keep5, "grid")
expect_equal(dim(grid5), c(18, 4))
expect_equal(grid5$thr.imp, rep(1, 18))
expect_equal(grid5$grp.num, rep(1:3, each = 6))

################################################################################

expect_error(snp_grid_PRS(G, all_keep, betas, lpval, type = "integer"))

multi_PRS <- snp_grid_PRS(G, all_keep, betas, lpval, type = "double",
                          n_thr_lpS = (n <- sample(10:30, 1)))
expect_identical(typeof(multi_PRS), "double")
expect_equal(dim(multi_PRS), c(nrow(G), n * sum(lengths(all_keep))))

multi_PRS <- snp_grid_PRS(G, all_keep, betas, lpval, grid.lpS.thr = 0:5, ncores = 2)
expect_identical(typeof(multi_PRS), "float")
expect_equal(dim(multi_PRS), c(nrow(G), 6 * sum(lengths(all_keep))))

multi_PRS2 <- lapply(all_keep2, function(ind.keep) {
  snp_PRS(G, betas[ind.keep], ind.keep = ind.keep, lpS.keep = lpval[ind.keep],
          thr.list = 0:5)
})
expect_equal(do.call("cbind", multi_PRS2),
             multi_PRS[, 1:36] + multi_PRS[, 37:72],
             check.attributes = FALSE, tolerance = 1e-7)

multi_PRS3 <- lapply(unlist(all_keep, recursive = FALSE), function(ind.keep) {
  snp_PRS(G, betas[ind.keep], ind.keep = ind.keep, lpS.keep = lpval[ind.keep],
          thr.list = 0:5)
})
expect_equal(do.call("cbind", multi_PRS3), multi_PRS[],
             check.attributes = FALSE, tolerance = 1e-7)

################################################################################

new_betas <- snp_grid_stacking(multi_PRS, y, alphas = 1e-3, ncores = 2)
expect_length(new_betas$beta.covar, 0)
expect_equal(
  predict(new_betas$mod, multi_PRS, proba = FALSE),
  new_betas$intercept + big_prodVec(G, new_betas$beta.G),
  check.attributes = FALSE, tolerance = 1e-6)

################################################################################
