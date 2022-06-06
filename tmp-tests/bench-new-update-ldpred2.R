library(bigsnpr)
options(max.print = 500)

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

corr <- runonce::save_run({
  POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data")
  plot(POS, POS2, pch = 20)
  snp_cor(chr6$genotypes, infos.pos = POS2, size = 3 / 1000, ncores = 6)
}, file = "tmp-data/corr_chr6.rds")

corr2 <- bigsparser::as_SFBM(corr, compact = TRUE)

# Simu phenotype
set.seed(1)
simu <- snp_simuPheno(G, h2 = 0.1, M = 100)
simu$effects
y2 <- simu$pheno


# GWAS
ind.gwas <- sample(nrow(G), 8e3)
gwas <- big_univLinReg(G, y2[ind.gwas], ind.train = ind.gwas)
# plot(gwas, type = "Manhattan")

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDSc reg
df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err, n_eff = length(ind.gwas))
(ldsc <- snp_ldsc2(corr, df_beta))

(h2_seq <- round(ldsc[["h2"]] * c(0.3, 1.4), 4))
(p_seq <- signif(seq_log(1e-5, 1, length.out = 11), 2))
(params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))


beta_grid <- sapply(rows_along(params), function(ic) {
  print(ic)
  params$time[ic] <<- system.time(
    res <- snp_ldpred2_grid(corr2, df_beta, grid_param = params[ic, ])
  )[3]
  params$conv[ic] <<- !anyNA(res)
  res
})


ord <- with(params, order(-p, sparse, -h2))
plot(params$time, pch = 20)
plot(params$time[ord], pch = 20)
subset(params, p == 1)



pred_grid <- big_prodMat(G, beta_grid, ind.row = ind.val)
summary(drop(cor(pred_grid, y2[ind.val])))
## Before:
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 0.02348 0.08969 0.12714 0.11332 0.14006 0.15522       5
## After:
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 0.02311 0.08054 0.12710 0.10799 0.13632 0.15519       2

library(ggplot2)
ggplot(params, aes(x = p, y = ifelse(conv, time, NA), color = as.factor(h2))) +
  theme_bigstatsr() +
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  geom_hline(yintercept = 0) +
  facet_wrap(~ sparse, labeller = label_both) +
  labs(y = "Runtime", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))

sum(params$time) / (nrow(params) * 1.5)


system.time(
  auto <- snp_ldpred2_auto(corr2, df_beta, h2_init = ldsc[["h2"]], verbose = TRUE)
)
# 9 sec -> 2-5 sec (more impressive for larger data)
