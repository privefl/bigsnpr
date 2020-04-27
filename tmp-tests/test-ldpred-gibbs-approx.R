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
h2 <- 0.2; M <- 1000
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
m <- ncol(G)
coeff <- N * h2 / m
corr2 <- corr
diag(corr2) <- 1 + 1 / coeff

chi2 <- qchisq(predict(gwas) * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)
betas_hat <- sqrt(chi2) * sign(beta_gwas) / sqrt(N)
new_beta <- as.vector(Matrix::solve(corr2, betas_hat))
pred <- big_prodVec(G, sqrt(N) * gwas$std.err * new_beta, ind.row = ind.val)
plot(pred, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred, y2[ind.val])**2


# LDpred-gibbs

curr_betas <- avg_betas <- rep(0, m)
p <- 0.1
post_p <- rep(p, m)
h2 <- 1
h2_est <- c()
burn_in <- 5000
# ord <- order(new_beta ** 2, decreasing = TRUE)

for (k in seq_len(50)) {

  coeff <- N * h2 / m
  diag(corr2) <- 1 + 1 / coeff
  new_beta <- as.vector(Matrix::solve(corr2, betas_hat))

  C <- coeff / (coeff + 1)
  Vnc <- C / N * (1 - C * h2)
  C2 <- (coeff + 1) / (coeff + p)
  Vc <- Vnc + C * h2 / (m * p)

  post_p <- 1 / (1 + (1 - p) / p * sqrt(Vc / Vnc) *
                   exp(0.5 * new_beta ** 2 * (1 / Vc - 1 / Vnc)))
  curr_betas <-
    ifelse(post_p > runif(m), C2 * new_beta + rnorm(m, sd = sqrt(C2 / N)), 0)

  if (k %% 2) {
    h2_est[k] <- drop(crossprod(curr_betas))
    h2 <- min(max(1e-4, (h2^3 * h2_est[k])^(1/4)), 1)
  } else {
    p <- max(1e-4, (p^3 * mean(post_p))^(1/4))
  }
  print(c(h2, p))

  # if (k > burn_in) avg_betas <- avg_betas + curr_post_means
}



avg_betas <- avg_betas / (num_iter - burn_in)
plot(avg_betas, new_beta, pch = 20); abline(0, 1, col = "red")

avg_betas2 <-  C2 * new_beta /
  (1 + (1 - p) / p * sqrt(Vc / Vnc) *
     exp(0.5 * new_beta ** 2 * (1 / Vc - 1 / Vnc)))
plot(avg_betas, avg_betas2); abline(0, 1, col = "red")

pred5 <- big_prodVec(G, sqrt(N) * gwas$std.err * avg_betas, ind.row = ind.val)
plot(pred5, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred5, y2[ind.val])**2
cor(pred, y2[ind.val])**2

plot(h2_est); abline(h = true_h2, col = "red")

