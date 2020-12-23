################################################################################

context("PRS")

################################################################################

test <- snp_attachExtdata()
G <- test$genotypes
y01 <- test$fam$affection - 1

################################################################################

# test_that()

# PCA -> covariables
obj.svd <- snp_autoSVD(G, infos.chr = test$map$chromosome,
                       infos.pos = test$map$physical.pos, verbose = FALSE)

# GWAS
gwas <- big_univLogReg(G, y01.train = y01, covar.train = obj.svd$u)
pval <- predict(gwas, log10 = FALSE)
pval2 <- readRDS(system.file("testdata", "pval.rds", package = "bigsnpr"))
expect_equal(pval, pval2, tolerance = 1e-4)

# clumping
ind.keep <- snp_clumping(G, infos.chr = test$map$chromosome,
                         S = abs(gwas$score),
                         size = 250, # as PLINK default
                         infos.pos = test$map$physical.pos)
ind.keep2 <- readRDS(system.file("testdata", "clumping.rds",
                                 package = "bigsnpr"))
expect_gt(mean(ind.keep %in% ind.keep2), 0.98)

# PRS
thrs <- seq(0, 5, by = 0.5)
prs <- snp_PRS(G, betas.keep = gwas$estim[ind.keep],
               ind.test = rows_along(G),
               ind.keep = ind.keep,
               lpS.keep = -predict(gwas)[ind.keep],
               thr.list = thrs)
expect_equal(dim(prs), c(nrow(G), length(thrs)))
prs2 <- readRDS(system.file("testdata", "scores-PRS.rds", package = "bigsnpr"))
scores.cor <- sapply(cols_along(prs), function(j) cor(prs[, j], prs2[, j]))
expect_equal(scores.cor, rep(1, length(thrs)),
             tolerance = 1e-3)

# No ordering in `thrs`
thrs2 <- sample(thrs)
prs <- snp_PRS(G, betas.keep = gwas$estim[ind.keep],
               ind.test = rows_along(G),
               ind.keep = ind.keep,
               lpS.keep = -predict(gwas)[ind.keep],
               thr.list = thrs2)
prs3 <- prs[, order(thrs2)]
expect_equal(dim(prs3), c(nrow(G), length(thrs2)))
scores.cor2 <- sapply(cols_along(prs3), function(j) cor(prs3[, j], prs2[, j]))
expect_equal(scores.cor2, rep(1, length(thrs2)),
             tolerance = 1e-3)

# No threshold
expect_message(snp_PRS(G, betas.keep = gwas$estim[ind.keep],
                       ind.test = rows_along(G),
                       ind.keep = ind.keep,
                       lpS.keep = -predict(gwas)[ind.keep]),
               "Thresholding disabled.")
expect_message(snp_PRS(G, betas.keep = gwas$estim[ind.keep],
                       ind.test = rows_along(G),
                       ind.keep = ind.keep,
                       thr.list = thrs2),
               "Thresholding disabled.")

################################################################################

test_that("snp_thr_correct() works", {

  beta <- rnorm(1000)
  beta_se <- runif(1000, min = 0.3, max = 0.5)
  lpval <- -log10(pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE))

  THR <- runif(1, 0, 2)
  new_beta <- snp_thr_correct(beta, beta_se = beta_se, thr_lpS = THR)
  new_beta2 <- snp_thr_correct(beta, lpS = lpval, thr_lpS = THR)
  expect_equal(new_beta2, new_beta)
  expect_error(snp_thr_correct(beta, thr_lpS = THR), "cannot be both missing")

  is_signif <- (lpval >= THR)
  expect_true(all(new_beta2[is_signif] != 0))
  expect_true(all(new_beta2[!is_signif] == 0))

  is_high_signif <- (lpval > 10)
  expect_equal(new_beta2[is_high_signif], beta[is_high_signif], tolerance = 1e-4)

  expect_gt(cor((new_beta2 / beta_se)[is_signif], (beta / beta_se)[is_signif],
                method = "spearman"), 0.999)
  expect_gt(cor((new_beta / beta)[is_signif], abs(beta / beta_se)[is_signif],
                method = "spearman"), 0.999)
  expect_true(all(sign(new_beta[is_signif]) == sign(beta[is_signif])))
  expect_true(all(abs(new_beta / beta_se) <= abs(beta / beta_se)))
})


################################################################################
