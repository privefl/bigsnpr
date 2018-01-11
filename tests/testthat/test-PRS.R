################################################################################

context("PRS")

################################################################################

test <- snp_attachExtdata()
G <- test$genotypes
y01 <- test$fam$affection - 1

################################################################################

# PCA -> covariables
obj.svd <- snp_autoSVD(G, infos.chr = test$map$chromosome,
                       infos.pos = test$map$physical.position)

# GWAS
gwas <- big_univLogReg(G, y01.train = y01, covar.train = obj.svd$u)
pval <- predict(gwas, log10 = FALSE)
pval2 <- readRDS(system.file("testdata", "pval.rds", package = "bigsnpr"))
expect_equal(pval, pval2, tolerance = 1e-4)

# clumping
ind.keep <- snp_clumping(G, infos.chr = test$map$chromosome,
                         S = abs(gwas$score),
                         size = 250, # as PLINK default
                         is.size.in.bp = TRUE,
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

################################################################################
