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

map_chr6 <- bigreadr::fread2("https://github.com/joepickrell/1000-genomes-genetic-maps/raw/master/interpolated_OMNI/chr6.OMNI.interpolated_genetic_map.gz")
ind <- bigutilsr::knn_parallel(as.matrix(map_chr6$V2), as.matrix(POS),
                               k = 1, ncores = 1)$nn.idx
POS2 <- map_chr6$V3[ind]
plot(POS, POS2, pch = 20)

corr <- snp_cor(chr6$genotypes, infos.pos = POS2, size = 3 / 1000,
                ncores = 6, alpha = 1)
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
h2 <- 0.05
M <- 2000; set <- sort(sample(ncol(G), size = M))
# M <- 10; set <- sort(sample(ind.HLA, size = M))
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
beta_est <- as.vector(Matrix::solve(corr2, betas_hat))
new_beta <- sqrt(N) * gwas$std.err * beta_est
pred <- big_prodVec(G, new_beta, ind.row = ind.val)
plot(pred, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred, y2[ind.val])**2  # 0.214 (2Mb) -> 0.256 (5Mb) -> 0.259 (5cM)

# LDpred-gibbs

snp_ldsc(ld / ncol(corr), chi2, N, blocks = 20)

Rcpp::sourceCpp('src/ldpred2-auto.cpp')
# burn_in <- 2000
nrep <- 1
all_ldpred <- ldpred2_gibbs_auto(
  corr      = corr,
  betas_hat = betas_hat,
  n_vec     = `if`(length(N) == 1, rep(N, m), N),
  h2_init   = rep(0.1, nrep),
  p_init    = seq_log(1e-4, 0.9, length.out = nrep),
  burn_in   = 2000,
  num_iter  = 1000,
  ncores    = min(6, nrep)
)
c(M / ncol(G), var(y) / var(y2))


ldpred <- all_ldpred[[1]]
str(ldpred)

# plot(ldpred$beta_est, betas_init, pch = 20); abline(0, 1, col = "red")

new_beta2 <- sqrt(N) * gwas$std.err * ldpred$beta_est
plot(beta_gwas, new_beta2, pch = 20); abline(0, 1, col = "red")
# plot(new_beta, new_beta2, pch = 20)
pred7 <- big_prodVec(G, new_beta2, ind.row = ind.val)
plot(pred7, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
round(100 * drop(cor(pred7, y2[ind.val])**2), 2)
round(100 * drop(cor(pred, y2[ind.val])**2), 2)

plot(ldpred$vec_p_est,  log = "y", pch = 20)
abline(h = print(ldpred$p_est),  col = "red", lwd = 2)
abline(h = print(M / ncol(G)), col = "blue", lwd = 2)
p_est <- tail(ldpred$vec_p_est, -1000)
abline(h = quantile(p_est, probs = c(0.025, 0.975)), col = "green", lwd = 2)

plot(ldpred$vec_h2_est, log = "y", pch = 20)
abline(h = print(ldpred$h2_est), col = "red", lwd = 2)
abline(h = print(var(y) / var(y2)), col = "blue", lwd = 2)
h2_est <- tail(ldpred$vec_h2_est, -1000)
abline(h = quantile(h2_est, probs = c(0.025, 0.975)), col = "green", lwd = 2)
snp_ldsc(ld / ncol(corr), chi2, N, blocks = 50)

