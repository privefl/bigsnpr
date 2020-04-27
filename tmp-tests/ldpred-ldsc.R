library(bigsnpr)

chr22 <- snp_attach("../Dubois2010_data/celiac_chr22.rds")
G <- chr22$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
dim(G)  # 11402  4945
big_counts(G, ind.col = 1:10)
CHR <- chr22$map$chromosome
POS <- chr22$map$physical.pos

corr <- snp_cor(chr22$genotypes, infos.pos = POS, size = 3000,
                ncores = 4, alpha = 1)
# ld <- Matrix::colMeans(corr ** 2 - (1 - corr ** 2) / (nrow(G) - 2))
ld <- Matrix::colMeans(corr ** 2)
plot(ld, pch = 20)
plot(bigutilsr::rollmean(ld, 100), pch = 20)

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
N <- length(ind.gwas)
gwas <- big_univLinReg(G, y2[ind.gwas], ind.train = ind.gwas)
beta_gwas <- gwas$estim

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDSC
chi2 <- qchisq(predict(gwas) * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)
(ldsc <- snp_ldsc(ld_div_size = ld, chi2 = chi2, sample_size = N,
                  chi2_thr1 = 30))
h2_ldsc <- min(max(0.01, ldsc[3]), 1)

sapply(setNames(nm = 2:8 * 10), function(thr) {
  snp_ldsc(ld_div_size = ld, chi2 = chi2, sample_size = N, chi2_thr1 = thr,
           blocks = 5)
})

# LDpred-inf
m <- length(beta_gwas)
coeff <- N * h2_ldsc / m
corr2 <- corr + Matrix::Diagonal(ncol(corr), 1 / coeff)
new_beta <- as.vector(Matrix::solve(corr2, beta_gwas))
pred <- big_prodVec(G, new_beta, ind.row = ind.val)
plot(pred, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred, y2[ind.val])**2  # 0.376 -> 0.38


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

curr_betas <- beta_hats <- sqrt(chi2) * sign(beta_gwas) / sqrt(N)
curr_post_means <- avg_betas <- rep(0, m)

burn_in <- 20
num_iter <- 100

for (k in seq_len(num_iter)) {

  # print(k)
  print(h2_est <- max(0.00001, crossprod(curr_betas)))
  alpha <- 1 #min(0.99, 1 / h2_est, (h2 + 1 / sqrt(N)) / h2_est)

  for (i in seq_len(m)) {
    curr_betas[i] <- 0
    res_beta_hat_i <- beta_hats[i] -
      crossprod(corr_as_list[[i]][[2]], curr_betas[corr_as_list[[i]][[1]]])
    postp <- 1 / (1 + C1 * exp(C3 * res_beta_hat_i ** 2))
    curr_post_means[i] <- C2 * postp * res_beta_hat_i
    curr_betas[i] <- if (runif(1) < (postp * alpha)) {
      rnorm(1) * C4 + C2 * res_beta_hat_i
    } else 0
  }

  if (k > burn_in) avg_betas <- avg_betas + curr_post_means
}

avg_betas <- avg_betas / (num_iter - burn_in)
plot(avg_betas, beta_hats, pch = 20); abline(0, 1, col = "red", lwd = 3)
pred4 <- big_prodVec(G, avg_betas, ind.row = ind.val)
plot(pred4, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred4, y2[ind.val])**2  # 0.348 -> 0.40
# completely off when
# - window size is too small (500 Kb or even 1 Mb, but okay with 2 Mb)
# - introduce strong sparsity (alpha = 1e-10, but okay with 0.01)

ldpred_shrink <- avg_betas / beta_hats
pred5 <- big_prodVec(G, beta_gwas * ldpred_shrink, ind.row = ind.val)
plot(pred5, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred5, y2[ind.val])**2
