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
# snp_subset(celiac, ind.row = which(log(dist) < 4), ind.col = which(CHR == 22),
#            backingfile = "../Dubois2010_data/celiac_chr22")

chr22 <- snp_attach("../Dubois2010_data/celiac_chr22.rds")
G <- chr22$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
dim(G)  # 11402  4945
big_counts(G, ind.col = 1:10)
CHR <- chr22$map$chromosome
POS <- chr22$map$physical.pos


# download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/1kg_eur.tar.bz2",
#               destfile = "tmp-data/1kg_eur.tar.bz2")
# genetic_map <- bigreadr::fread2("tmp-data/22.bim")
# ind <- bigutilsr::knn_parallel(as.matrix(genetic_map$V4), as.matrix(POS), k = 1)$nn.idx
# POS2 <- genetic_map$V3[ind]

corr <- snp_cor(chr22$genotypes, infos.pos = POS, size = 5000,
                ncores = 4, alpha = 1)
str(corr)
median(Matrix::colSums(corr != 0))

# Simu phenotype
set.seed(1)
h2 <- 0.5; M <- 1000
set <- sort(sample(ncol(G), size = M))
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
coeff <- N * h2 / ncol(G)
corr2 <- corr + Matrix::Diagonal(ncol(corr), 1 / coeff)
new_beta <- as.vector(Matrix::solve(corr2, beta_gwas))
pred <- big_prodVec(G, new_beta, ind.row = ind.val)
plot(pred, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred, y2[ind.val])**2  # 0.376 -> 0.38

# Naive PRS
pred2 <- big_prodVec(G, beta_gwas, ind.row = ind.val)
plot(pred2, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred2, y2[ind.val])**2  # 0.288

# PRS with clumping
ind.keep <- snp_clumping(G, infos.chr = CHR, infos.pos = POS, ind.row = ind.val,
                         S = -predict(gwas), ncores = 4)
pred3 <- big_prodVec(G, beta_gwas[ind.keep], ind.row = ind.val, ind.col = ind.keep)
plot(pred3, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred3, y2[ind.val])**2  # 0.318

# LDpred-gibbs
corr2 <- as(corr, "dgTMatrix")
corr_as_list <- split(data.frame(i = corr2@i + 1L, r = corr2@x),
                      factor(corr2@j, ordered = TRUE))

p <- 0.2
L <- coeff / p
C1 <- (1 - p) / p * sqrt(1 + L)
C2 <- L / (L + 1)
C3 <- -N / 2 * C2
C4 <- sqrt(C2 / N)

m <- length(beta_gwas)
chi2 <- qchisq(predict(gwas) * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)
beta_hats <- sqrt(chi2) * sign(beta_gwas) / sqrt(N)
curr_betas <- beta_hats / 1   # should start with beta_gwas instead?
curr_post_means <- rep(0, m)
avg_betas <- rep(0, m)

burn_in <- 20
num_iter <- 100

all_h2_est <- c()

for (k in seq_len(num_iter)) {

  # print(k)
  print(all_h2_est[k] <- h2_est <- max(0.00001, crossprod(curr_betas)))
  print(crossprod(curr_betas, corr %*% curr_betas))
  alpha <- 1 #min(0.99, 1 / h2_est, (h2 + 1 / sqrt(N)) / h2_est)

  for (i in seq_len(m)) {
    curr_betas[i] <- 0
    res_beta_hat_i <- beta_hats[i] -
      crossprod(corr_as_list[[i]][[2]], curr_betas[corr_as_list[[i]][[1]]])
    postp <- 1 / (1 + C1 * exp(C3 * res_beta_hat_i ** 2))
    curr_post_means[i] <- C2 * postp * res_beta_hat_i
    curr_betas[i] <- if (runif(1) < (postp * alpha)) {
      rnorm(1, sd = C4) + C2 * res_beta_hat_i
    } else 0
  }

  if (k > burn_in) avg_betas <- avg_betas + curr_post_means
}

avg_betas <- avg_betas / (num_iter - burn_in)
plot(avg_betas, new_beta, pch = 20)
pred4 <- big_prodVec(G, avg_betas, ind.row = ind.val)
plot(pred4, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred4, y2[ind.val])**2  # 0.348 -> 0.40
# completely off when
# - window size is too small (500 Kb or even 1 Mb, but okay with 2 Mb)
# - introduce strong sparsity (alpha = 1e-10, but okay with 0.01)

ldpred_shrink <- avg_betas / beta_hats
# shrink_max <- 5
# ldpred_shrink2 <- pmin(pmax(-shrink_max, ldpred_shrink), shrink_max)
pred5 <- big_prodVec(G, beta_gwas * ldpred_shrink, ind.row = ind.val)
plot(pred5, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred5, y2[ind.val])**2

plot(all_h2_est)
