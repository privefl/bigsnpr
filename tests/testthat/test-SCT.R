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
  check.attributes = FALSE, tolerance = 1e-7)

################################################################################

sumstats <- data.frame(
  chr = 1,
  pos = c(86303, 86331, 162463, 752566, 755890, 758144),
  a0 = c("T", "G", "C", "A", "T", "G"),
  a1 = c("G", "A", "T", "G", "A", "A"),
  beta = c(-1.868, 0.250, -0.671, 2.112, 0.239, 1.272),
  p = c(0.860, 0.346, 0.900, 0.456, 0.776, 0.383)
)

info_snp <- data.frame(
  id = c("rs2949417", "rs115209712", "rs143399298", "rs3094315", "rs3115858"),
  chr = 1,
  pos = c(86303, 86331, 162463, 752566, 755890),
  a0 = c("T", "A", "G", "A", "T"),
  a1 = c("G", "G", "A", "G", "A")
)

expect_message(matched1 <- snp_match(sumstats, info_snp),
               "4 variants have been matched; 1 were flipped and 1 were reversed.")
expect_equal(dim(matched1), c(4, 7))
expect_equal(matched1$beta, sumstats$beta[1:4] * c(1, -1, 1, 1))
expect_message(matched2 <- snp_match(sumstats, info_snp, strand_flip = FALSE),
               "4 variants have been matched; 0 were flipped and 1 were reversed.")
expect_equal(dim(matched2), c(4, 7))
expect_equal(matched2$beta, sumstats$beta[c(1:2, 4:5)] * c(1, -1, 1, 1))

################################################################################
