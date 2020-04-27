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
corr3 <- as(corr, "dgTMatrix")
corr_as_list <- split(data.frame(i = corr3@i + 1L, r = corr3@x),
                      factor(corr3@j, ordered = TRUE))
object.size(corr_as_list) / 1024**2 # 372 Mb

chi2 <- qchisq(predict(gwas) * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)
beta_hats <- sqrt(chi2) * sign(beta_gwas) / sqrt(N)

sapply(setNames(nm = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3)), function(p) {

  print(p)

  L <- coeff / p
  C1 <- (1 - p) / p * sqrt(1 + L)
  C2 <- L / (L + 1)
  C3 <- -N / 2 * C2
  C4 <- sqrt(C2 / N)

  curr_betas <- curr_post_means <- avg_betas <- rep(0, m)
  burn_in <- 10
  num_iter <- 60

  for (k in seq_len(num_iter)) {

    # print(k)
    # print(h2_est <- max(0.00001, crossprod(curr_betas)))
    alpha <- 1 #min(0.99, 1 / h2_est, (h2 + 1 / sqrt(N)) / h2_est)

    for (i in seq_len(m)) {
      curr_betas[i] <- 0
      res_beta_hat_i <- beta_hats[i] -
        crossprod(corr_as_list[[i]][[2]], curr_betas[corr_as_list[[i]][[1]]])
      postp <- 1 / (1 + C1 * exp(C3 * res_beta_hat_i ** 2))
      if (postp < p) {
        curr_betas[i] <- curr_post_means[i] <- 0
      } else {
        curr_post_means[i] <- C2 * postp * res_beta_hat_i
        curr_betas[i] <- if (runif(1) < (postp * alpha)) {
          rnorm(1) * C4 + C2 * res_beta_hat_i
        } else 0
      }
    }

    if (k > burn_in) avg_betas <- avg_betas + curr_post_means
  }

  avg_betas <- avg_betas / (num_iter - burn_in)
  ldpred_shrink <- avg_betas / beta_hats

  pred5 <- big_prodVec(G, beta_gwas * ldpred_shrink, ind.row = ind.val)
  round(100 * c(cor(pred5, y2[ind.val])**2, mean(avg_betas == 0)), 2)
})
## postp < p
#      0.001 0.003  0.01  0.03   0.1   0.3
# [1,] 17.30 21.47 25.67 29.53 31.13 29.84
# [2,] 84.36 78.34 65.76 50.41 34.16 28.62
## nothing
#      0.001 0.003  0.01  0.03   0.1   0.3
# [1,] 17.82 21.89 26.12 30.27 31.18 28.62
# [2,]  0.00  0.00  0.00  0.00  0.00  0.00
## postp < (p / 2)
#      0.001 0.003  0.01  0.03   0.1   0.3
# [1,] 17.52 21.35 26.16 29.60 30.92 28.62
# [2,] 77.51 65.13 42.48 10.01  0.00  0.00
