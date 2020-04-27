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
h2 <- 0.1; M <- 100
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
chi2 <- qchisq(predict(gwas) * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDpred-inf
m <- length(beta_gwas)
coeff <- N * h2 / m
corr2 <- corr + Matrix::Diagonal(ncol(corr), 1 / coeff)
betas_blup <- as.vector(Matrix::solve(corr2, beta_gwas))
pred <- big_prodVec(G, betas_blup, ind.row = ind.val)
plot(pred, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred, y2[ind.val])**2  # 0.376 -> 0.38


# LDpred-gibbs
corr3 <- as(corr, "dgTMatrix")
corr_as_list <- split(data.frame(i = corr3@i + 1L, r = corr3@x),
                      factor(corr3@j, ordered = TRUE))

sapply(setNames(nm = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3)), function(p) {

  print(p)

  # fast as starting values
  C1 <- (coeff + 1) / (coeff + p)
  C2 <- coeff / (coeff + 1)
  V_nc <- C2 / N * (1 - C2 * h2)
  V_c <- C2 * h2 / (m * p) + V_nc
  d_beta_ncaus <- dnorm(betas_blup, sd = sqrt(V_nc))
  d_beta_caus <- dnorm(betas_blup, sd = sqrt(V_c))
  d_caus_beta <- d_beta_caus * p / (d_beta_caus * p + d_beta_ncaus * (1 - p))

  beta_hats <- sqrt(chi2) * sign(beta_gwas) / sqrt(N)
  betas_blup <- as.vector(Matrix::solve(corr2, beta_hats))
  fast_betas <- betas_blup * C1 * d_caus_beta
  curr_betas <- beta_hats #fast_betas

  L <- coeff / p
  C1 <- (1 - p) / p * sqrt(1 + L)
  C2 <- L / (L + 1)
  C3 <- -N / 2 * C2
  C4 <- sqrt(C2 / N)

  curr_post_means <- avg_betas <- rep(0, m)
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
      curr_post_means[i] <- C2 * postp * res_beta_hat_i
      curr_betas[i] <- if (runif(1) < (postp * alpha)) {
        rnorm(1) * C4 + C2 * res_beta_hat_i
      } else 0
    }

    if (k > burn_in) avg_betas <- avg_betas + curr_post_means
  }

  avg_betas <- avg_betas / (num_iter - burn_in)
  ldpred_shrink <- avg_betas / beta_hats

  plot(avg_betas, fast_betas, pch = 20)

  pred5 <- big_prodVec(G, beta_gwas * ldpred_shrink, ind.row = ind.val)
  # plot(pred5, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
  round(100 * cor(pred5, y2[ind.val])**2, 1)
})

# Naive PRS
pred2 <- big_prodVec(G, beta_gwas, ind.row = ind.val)
plot(pred2, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred2, y2[ind.val])**2

# PRS with clumping
ind.keep <- snp_clumping(G, infos.chr = CHR, infos.pos = POS, ind.row = ind.val,
                         S = -predict(gwas), ncores = 4)
pred3 <- big_prodVec(G, beta_gwas[ind.keep], ind.row = ind.val, ind.col = ind.keep)
plot(pred3, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred3, y2[ind.val])**2
