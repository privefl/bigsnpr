library(bigsnpr)

# celiac <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")
# G <- celiac$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
# CHR <- celiac$map$chromosome
# POS <- celiac$map$physical.pos
#
# obj.svd <- snp_autoSVD(G, CHR, POS, ncores = 4, thr.r2 = 0.1)
# plot(obj.svd)
#
# dist <- bigutilsr::dist_ogk(obj.svd$u)
# hist(log(dist))
#
# snp_subset(celiac, ind.row = which(log(dist) < 4), ind.col = which(CHR == 6),
#            backingfile = "../Dubois2010_data/celiac_chr6")

chr6 <- snp_attach("../Dubois2010_data/celiac_chr6.rds")
G <- chr6$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
dim(G)
big_counts(G, ind.col = 1:10)
CHR <- chr6$map$chromosome
POS <- chr6$map$physical.pos

POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data/")
plot(POS, POS2, pch = 20)

corr <- runonce::save_run(
  snp_cor(chr6$genotypes, infos.pos = POS2, size = 3 / 1000, ncores = 6),
  file = "tmp-data/corr_chr6.rds"
)

# Simu phenotype
ind.HLA <- snp_indLRLDR(CHR, POS, LD.wiki34[12, ])
set.seed(1)
# simu <- snp_simuPheno(G, h2 = 0.5, M = 1000, ind.possible = ind.HLA)
# y2 <- simu$pheno
y2 <- snp_simuPheno(G, h2 = 0.2, M = 500)$pheno
g <- big_prodVec(G, rnorm(500, sd = sqrt(0.2 / 500)),
                 ind.col = sample(ncol(G), 500))
y2 <- g / sd(g) * sqrt(0.2) + rnorm(length(g), sd = sqrt(0.8))


# GWAS
ind.gwas <- sample(nrow(G), 8e3)
gwas <- big_univLinReg(G, y2[ind.gwas], ind.train = ind.gwas)
plot(gwas, type = "Manhattan")

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDSc reg
df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err, n_eff = length(ind.gwas))
(ldsc <- snp_ldsc2(corr, df_beta))
h2_est <- ldsc[["h2"]]

corr2 <- bigsparser::as_SFBM(corr)

# LDpred-auto
auto <- snp_ldpred2_auto(corr2, df_beta, h2_init = h2_est,
                         burn_in = 200, num_iter = 200,
                         report_step = 5,
                         verbose = TRUE)[[1]]
# Overall: 0.0206819 // 0.200177

median(alpha_est <- tail(auto$path_alpha_est, 200) - 1)
hist(alpha_est)

pred <- big_prodVec(G, auto$beta_est, ind.row = ind.val)
cor(pred, y2[ind.val])

bsamp <- auto$sample_beta
sd <- with(df_beta, 1 / sqrt(n_eff * beta_se^2 + beta^2))
y <- bsamp[, 2]
ind <- which(y != 0)
plot(log(y[ind]^2) ~ log(sd[ind]))
summary(lm(log(y[ind]^2) ~ log(sd[ind])))$coef[2, 1:2]
alpha_est[1:10] + 1


all_coef <- apply(bsamp, 2, function(y) {
  ind <- which(y != 0)
  summary(lm(log(y[ind]^2) ~ log(sd[ind])))$coef[2, 1:2]
})
hist(all_coef[1, ] - 1)

all_coef2 <- apply(bsamp, 2, function(y) {
  ind <- which(y != 0)
  summary(MASS::rlm(log(y[ind]^2) ~ log(sd[ind])))$coef[2, 1:2]
})
hist(all_coef[1, ] - 1)
hist(all_coef2[1, ] - 1)
plot(all_coef[1, ], all_coef2[1, ])
