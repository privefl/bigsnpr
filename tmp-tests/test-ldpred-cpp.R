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
ind <- bigutilsr::knn_parallel(as.matrix(map_chr6$V2), as.matrix(POS), k = 1)$nn.idx
POS2 <- map_chr6$V3[ind]
plot(POS, POS2, pch = 20)

corr <- snp_cor(chr6$genotypes, infos.pos = POS2, size = 5 / 1000,
                ncores = 4, alpha = 1)
str(corr)
median(Matrix::colSums(corr != 0))
object.size(corr) / 1024**2

ld <- Matrix::colSums(corr ** 2)
plot(ld, pch = 20)
plot(POS, bigutilsr::rollmean(ld, 100))
hist(S <- bigutilsr::rollmean(ld, 100), "FD")
abline(v = (q <- bigutilsr::tukey_mc_up(S)), col = "red")

# Simu phenotype
ind.HLA <- snp_indLRLDR(CHR, POS, LD.wiki34[12, ])
set.seed(1)
h2 <- 0.5
M <- 1000; set <- sort(sample(ncol(G), size = M))
M <- 100; set <- sort(sample(ind.HLA, size = M))
effects <- rnorm(M, sd = sqrt(h2 / M))
y <- drop(scale(G[, set]) %*% effects)       ## G
y2 <- y + rnorm(nrow(G), sd = sqrt(1 - h2))  ## G + E
var(y) / var(y2)                             ## H2

# GWAS
ind.gwas <- sample(nrow(G), 8e3)
gwas <- big_univLinReg(G, y2[ind.gwas], ind.train = ind.gwas)
beta_gwas <- gwas$estim
chi2 <- qchisq(predict(gwas) * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)

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
p_vec <- seq_log(from = 0.001, to = 1, length.out = 20)

Rcpp::sourceCpp('tmp-tests/test-ldpred-cpp.cpp')
all_shrink <- ldpred_gibbs(
  corr      = corr,
  betas_hat = sqrt(chi2) * sign(beta_gwas) / sqrt(N),
  n_vec     = `if`(length(N) == 1, rep(N, m), N),
  p_vec     = p_vec,
  coeff     = coeff,
  burn_in   = 20,
  num_iter  = 100,
  sparse    = FALSE
)
# cor(all_shrink)

## For 6 values of p,
# R version: 270 sec and 68 GB temporary data
# Rcpp version: 23 sec

pred5 <- big_prodMat(G, sweep(all_shrink, 1, beta_gwas, '*'), ind.row = ind.val)
rbind(p = round(p_vec, 4),
      r2 = round(100 * drop(cor(pred5, y2[ind.val])**2), 2),
      sp = round(100 * colMeans(all_shrink == 0), 2))

#      [,1]    [,2]    [,3]   [,4]    [,5]    [,6]    [,7]    [,8]    [,9]   [,10]
# p   0.001  0.0014  0.0021  0.003  0.0043  0.0062  0.0089  0.0127  0.0183  0.0264
# r2 46.850 47.6000 48.9200 48.090 48.3500 48.2500 48.7100 48.6000 48.1200 47.6000
# sp  0.000  0.0000  0.0000  0.000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
#      [,11]   [,12]   [,13]   [,14]   [,15]   [,16]  [,17]   [,18]   [,19] [,20]
# p   0.0379  0.0546  0.0785  0.1129  0.1624  0.2336  0.336  0.4833  0.6952  1.00
# r2 46.1700 45.0000 43.5100 41.9800 40.2600 38.6600 37.080 35.9500 34.4100 33.23
# sp  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.000  0.0000  0.0000  0.00

#      [,1]    [,2]    [,3]   [,4]    [,5]    [,6]    [,7]    [,8]    [,9]   [,10]
# p   0.001  0.0014  0.0021  0.003  0.0043  0.0062  0.0089  0.0127  0.0183  0.0264
# r2 48.260 48.0500 48.2200 47.690 49.1800 48.7900 49.1000 48.5900 47.9300 47.4100
# sp 93.720 92.8800 92.3400 91.090 88.9800 86.6700 83.3700 80.0100 74.7100 69.5700
#      [,11]   [,12]   [,13]   [,14]   [,15]   [,16]  [,17]   [,18]   [,19] [,20]
# p   0.0379  0.0546  0.0785  0.1129  0.1624  0.2336  0.336  0.4833  0.6952  1.00
# r2 46.6100 45.9000 44.7200 43.7200 42.4000 41.4200 40.030 38.9400 37.5600 33.12
# sp 62.9900 56.5000 50.4000 45.8400 42.1500 39.9700 38.820 38.1900 37.6800  0.00


Rcpp::sourceCpp('tmp-tests/test-ldpred-cpp-postp.cpp')
shrink <- ldpred_gibbs_auto(
  corr      = corr,
  betas_hat = sqrt(chi2) * sign(beta_gwas) / sqrt(N),
  n_vec     = `if`(length(N) == 1, rep(N, m), N),
  p_init    = 0.1,
  coeff     = coeff,
  burn_in   = 100,
  num_iter  = 500,
  sparse    = FALSE
)

pred7 <- big_prodVec(G, beta_gwas * shrink, ind.row = ind.val)
round(100 * drop(cor(pred7, y2[ind.val])**2), 2)
mean(shrink == 0)
