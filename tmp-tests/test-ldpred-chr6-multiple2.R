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

corr <- snp_cor(chr6$genotypes, infos.pos = POS2, size = 3 / 1000,
                ncores = 6, alpha = 0.9)
object.size(corr) / 1024**2  # 0.3 -> 51 Mb / 0.9 -> 79 / 1 -> 84
str(corr)
median(Matrix::colSums(corr != 0))

ld <- Matrix::colSums(corr ** 2)
plot(ld, pch = 20)
plot(POS, bigutilsr::rollmean(ld, 100))
hist(S <- bigutilsr::rollmean(ld, 100), "FD")
abline(v = (q <- bigutilsr::tukey_mc_up(S)), col = "red")

# Simu phenotype
ind.HLA <- snp_indLRLDR(CHR, POS, LD.wiki34[12, ])
set.seed(1)
h2 <- 0.03
# M <- 500; set <- sort(sample(ncol(G), size = M))
M <- 1000; set <- sort(sample(ind.HLA, size = M))
# set <- ind.HLA[c(7, 8, 10, 12, 15)]; M <- 5
effects <- rnorm(M, sd = sqrt(h2 / M))
# effects <- rep(sqrt(h2 / M), M)
y <- drop(scale(G[, set]) %*% effects)       ## G
y2 <- y + rnorm(nrow(G), sd = sqrt(1 - h2))  ## G + E
var(y) / var(y2)                             ## H2


# GWAS
ind.gwas <- sample(nrow(G), 8e3)
gwas <- big_univLinReg(G, y2[ind.gwas], ind.train = ind.gwas)
beta_gwas <- gwas$estim
plot(gwas, type = "Manhattan")

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDpred-inf
N <- length(ind.gwas)
m <- ncol(G)
coeff <- N * h2 / m
corr2 <- corr + Matrix::Diagonal(ncol(corr), 1 / coeff)
sd <- big_scale()(G)$scale
chi2 <- qchisq(predict(gwas) * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)
betas_hat <- sqrt(chi2) * sign(beta_gwas) / sqrt(N)
betas_inf <- as.vector(Matrix::solve(corr2, betas_hat))
crossprod(betas_hat, betas_inf)

new_beta <- sqrt(N) * gwas$std.err * betas_inf
pred <- big_prodVec(G, new_beta, ind.row = ind.val)
plot(pred, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred, y2[ind.val])**2  # 0.214 (2Mb) -> 0.256 (5Mb) -> 0.259 (5cM)

# LDpred-gibbs

(ldsc <- snp_ldsc(ld, ncol(corr), chi2, N, blocks = 100))
(h2_seq <- round(ldsc[["h2"]] + -1:3 * ldsc[["h2_se"]], 3))
(p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2))
# 5 x 21 = 85

(params <- expand.grid(h2 = h2_seq, p = p_seq, sparse = c(FALSE)))


Rcpp::sourceCpp('src/ldpred2.cpp')
all_ldpred2 <- ldpred2_gibbs(
  corr       = corr,
  betas_hat  = betas_hat,
  betas_init = betas_inf,
  order      = order(betas_inf**2, decreasing = TRUE) - 1L,
  n_vec      = `if`(length(N) == 1, rep(N, m), N),
  h2         = params$h2,
  p          = params$p,
  sparse     = params$sparse,
  burn_in    = 100,
  num_iter   = 100,
  ncores     = 4
)

new_beta2 <- sweep(all_ldpred2, 1, sqrt(N) * gwas$std.err, '*')
pred7 <- big_prodMat(G, new_beta2, ind.row = ind.val)
params$r2 <- r2 <- drop(cor(pred7, y2[ind.val])**2)
round(100 * r2, 2)
round(100 * drop(cor(pred, y2[ind.val])**2), 2)
round(100 * drop(cor(pred7, pred)), 2)

library(ggplot2)
ggplot(params, aes(x = p, y = r2, color = factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  facet_wrap(~ sparse) +
  geom_vline(xintercept = M / ncol(G), linetype = 3)

library(dplyr)
params %>%
  mutate(sparsity = colMeans(all_ldpred2 == 0)) %>%
  arrange(desc(r2)) %>%
  mutate_at(4:5, round, digits = 3) %>%
  slice(1:10)
