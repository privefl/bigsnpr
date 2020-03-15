library(bigsnpr)

celiac <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")
G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos

base_url <- "https://github.com/joepickrell/1000-genomes-genetic-maps/raw/master/interpolated_OMNI/"
POS2 <- unlist(lapply(1:22, function(chr) {
  map_chr <- bigreadr::fread2(paste0(base_url, "chr", chr, ".OMNI.interpolated_genetic_map.gz"))
  ind_chr <- which(CHR == chr)
  ind <- bigutilsr::knn_parallel(as.matrix(map_chr$V2), as.matrix(POS[ind_chr]),
                                 k = 1, ncores = 1)$nn.idx
  map_chr$V3[ind]
}))

plot(POS, POS2,, col = CHR, pch = 20)

lrldr <- matrix(POS2[bigsnpr:::getIntervals(snp_indLRLDR(CHR, POS))],
                ncol = 2, byrow = F)
lrldr[, 2] - lrldr[, 1]

POS3 <- POS2 + 1000 * CHR

system.time(
  corr <- snp_cor(G, ind.row = sort(sample(nrow(G), 5e3)), size = 3 / 1000,
                  alpha = 1, infos.pos = POS3, ncores = 4)
)
#    user  system elapsed
# 1277.53  285.97  842.41


# Simu phenotype
G <- celiac$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
ind.HLA <- snp_indLRLDR(CHR, POS, LD.wiki34[12, ])
set.seed(1)
h2 <- 0.5
M <- 10e3; set <- sort(sample(ncol(G), size = M))
# M <- 100; set <- sort(sample(ind.HLA, size = M))
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

object.size(corr) / 1024**2  # 970 Mb
ld <- Matrix::colMeans(corr ** 2)

# LDpred-inf
N <- length(ind.gwas)
(ldsc <- snp_ldsc(ld_div_size = ld, chi2 = chi2, sample_size = N))
m <- ncol(G)
coeff <- N * ldsc[3] / m
corr2 <- corr + Matrix::Diagonal(ncol(corr), 1 / coeff)
new_beta <- as.vector(Matrix::solve(corr2, beta_gwas))
pred <- big_prodVec(G, new_beta, ind.row = ind.val)
plot(pred, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred, y2[ind.val])**2  # 3.2%

# LDpred-gibbs
p_vec <- seq_log(from = 0.0001, to = 1, length.out = 12)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
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

#     [,1]    [,2]    [,3]   [,4]    [,5]    [,6]    [,7]    [,8]   [,9]  [,10]
# p   0.001  0.0014  0.0021  0.003  0.0043  0.0062  0.0089  0.0127 0.0183 0.0264
# r2 35.860 29.7000 23.8800 18.800 15.5300 12.8800 11.2600 10.0600 9.2700 8.6600
# sp  0.000  0.0000  0.0000  0.000  0.0000  0.0000  0.0000  0.0000 0.0000 0.0000
#     [,11]  [,12]  [,13]  [,14]  [,15]  [,16] [,17]  [,18]  [,19] [,20]
# p  0.0379 0.0546 0.0785 0.1129 0.1624 0.2336 0.336 0.4833 0.6952  1.00
# r2 8.1700 7.7500 7.2200 6.6100 6.0700 5.4300 4.870 4.3000 3.7900  3.39
# sp 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.000 0.0000 0.0000  0.00

#      [,1]    [,2]    [,3]   [,4]    [,5]    [,6]    [,7]    [,8]    [,9]   [,10]
# p   0.001  0.0014  0.0021  0.003  0.0043  0.0062  0.0089  0.0127  0.0183  0.0264
# r2 36.430 30.8800 24.9200 20.420 16.4600 14.0900 12.3100 11.3100 10.6900 10.1400
# sp 85.440 81.9600 77.8500 73.140 67.9900 62.8200 58.0100 53.7300 50.4600 48.3200
#      [,11]   [,12]   [,13]   [,14]   [,15]   [,16]  [,17]   [,18]   [,19] [,20]
# p   0.0379  0.0546  0.0785  0.1129  0.1624  0.2336  0.336  0.4833  0.6952  1.00
# r2  9.7700  9.3200  8.9400  8.2400  7.6700  6.8800  6.170  5.5100  4.8600  3.39
# sp 46.8400 46.0800 45.6800 45.6900 45.6700 45.7700 45.980 46.1500 46.3200  0.00

Rcpp::sourceCpp('tmp-tests/test-ldpred-cpp-postp2.cpp')
shrink <- ldpred_gibbs_auto2(
  corr      = corr,
  betas_hat = sqrt(chi2) * sign(beta_gwas) / sqrt(N),
  n_vec     = `if`(length(N) == 1, rep(N, m), N),
  p_init    = 0.9,
  coeff     = coeff,
  burn_in   = 100,
  num_iter  = 500,
  sparse    = TRUE,
  w         = 1 / ld
) # Very slow to converge

pred7 <- big_prodVec(G, beta_gwas * shrink, ind.row = ind.val)
plot(pred7, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
round(100 * drop(cor(pred7, y2[ind.val])**2), 2)
mean(shrink == 0)
# Starting at 0.1:   0.422583 // 0.78387 // 4.58 // 0
# Starting at 0.001: 0.316331 // 0.789662 // 5.13 // 0
# With weights: 3.2695e-005 // 0.626209 // 53.3 // 97.4%
# With weights: 3.16633e-005 // 0.482992 // 54.14 // 0
