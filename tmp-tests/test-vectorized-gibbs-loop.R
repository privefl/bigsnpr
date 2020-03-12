library(bigsnpr)

chr22 <- snp_attach("../Dubois2010_data/celiac_chr22.rds")
G <- chr22$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
dim(G)  # 11402  4945
CHR <- chr22$map$chromosome
POS <- chr22$map$physical.pos

corr <- snp_cor(chr22$genotypes, infos.pos = POS, size = 2000,
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

# LDpred-gibbs
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
avg_betas <- rep(0, m)

burn_in <- 20
num_iter <- 100

postp <- rep(p, m)
true_beta <- rep(0, m); true_beta[set] <- effects

for (k in seq_len(num_iter)) {

  # print(k)
  print(h2_est <- max(0.00001, crossprod(curr_betas)))
  alpha <- 1 #min(0.99, 1 / h2_est, (h2 + 1 / sqrt(N)) / h2_est)

  corr2 <- corr + Matrix::Diagonal(ncol(corr), postp / coeff)
  res_beta_hat <- as.vector(Matrix::solve(corr2, curr_betas))
  postp <- 1 / (1 + C1 * exp(C3 * res_beta_hat ** 2))

  plot(postp, abs(true_beta))

  curr_post_means <- C2 * postp * res_beta_hat
  curr_betas <- ifelse(runif(m) < (postp * alpha),
                       rnorm(m, sd = C4) + C2 * res_beta_hat,
                       0)

  if (k > burn_in) avg_betas <- avg_betas + curr_post_means
}

avg_betas <- avg_betas / (num_iter - burn_in)
plot(avg_betas, new_beta, pch = 20)
pred4 <- big_prodVec(G, avg_betas, ind.row = ind.val)
plot(pred4, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred4, y2[ind.val])**2  # 0.348 -> 0.40

# DOES NOT WORK AT ALL LOL
