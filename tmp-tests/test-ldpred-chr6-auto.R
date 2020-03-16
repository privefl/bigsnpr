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
dim(G)  # 11402  4945
big_counts(G, ind.col = 1:10)
CHR <- chr6$map$chromosome
POS <- chr6$map$physical.pos

map_chr6 <- bigreadr::fread2("https://github.com/joepickrell/1000-genomes-genetic-maps/raw/master/interpolated_OMNI/chr6.OMNI.interpolated_genetic_map.gz")
ind <- bigutilsr::knn_parallel(as.matrix(map_chr6$V2), as.matrix(POS),
                               k = 1, ncores = 1)$nn.idx
POS2 <- map_chr6$V3[ind]
plot(POS, POS2, pch = 20)

corr <- snp_cor(chr6$genotypes, infos.pos = POS2, size = 3 / 1000,
                ncores = 4, alpha = 1)
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
h2 <- 0.1
M <- 100; set <- sort(sample(ncol(G), size = M))
M <- 1000; set <- sort(sample(ind.HLA, size = M))
effects <- rnorm(M, sd = sqrt(h2 / M))
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
new_beta <- as.vector(Matrix::solve(corr2, beta_gwas))
pred <- big_prodVec(G, new_beta, ind.row = ind.val)
plot(pred, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred, y2[ind.val])**2  # 0.214 (2Mb) -> 0.256 (5Mb) -> 0.259 (5cM)

# LDpred-gibbs
chi2 <- qchisq(predict(gwas) * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)

snp_ldsc(ld / ncol(corr), chi2, N, blocks = 50)

Rcpp::sourceCpp('tmp-tests/test-ldpred-cpp-postp3.cpp')
betas_hat <- sqrt(chi2) * sign(beta_gwas) / sqrt(N)
shrink <- ldpred_gibbs_auto3(
  corr      = corr,
  betas_hat = betas_hat,
  order     = order(betas_hat ** 2, decreasing = TRUE) - 1L,
  n_vec     = `if`(length(N) == 1, rep(N, m), N),
  h2_init   = 1,
  p_init    = 0.5,
  burn_in   = 200,
  num_iter  = 1000,
  sparse    = FALSE,
  w         = rep(1, length(ld))
) # h2_init = 1 works well // sparse = TRUE does not work when h2 can vary

pred7 <- big_prodVec(G, beta_gwas * shrink$shrink, ind.row = ind.val)
plot(pred7, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
round(100 * drop(cor(pred7, y2[ind.val])**2), 2)
round(100 * drop(cor(pred, y2[ind.val])**2), 2)

plot(shrink$vec_p_est,  log = "y", pch = 20)
abline(h = shrink$p_est,  col = "red", lwd = 2)
plot(shrink$vec_h2_est, log = "y", pch = 20)
abline(h = shrink$h2_est, col = "red", lwd = 2)

h2_est <- tail(shrink$vec_h2_est, 1000)
boot_h2_est <- replicate(10e3, {
  mean(sample(h2_est, replace = TRUE))
})
mean(h2_est)
# 0.4960173
quantile(boot_h2_est, probs = c(0.025, 0.975))
#      2.5%     97.5%
# 0.4947999 0.4972172
quantile(h2_est, probs = c(0.025, 0.975))
#      2.5%     97.5%
# 0.4580894 0.5360437
var(y) / var(y2)
