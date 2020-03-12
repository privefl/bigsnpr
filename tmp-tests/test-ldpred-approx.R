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

ind.val <- setdiff(rows_along(G), ind.gwas)

# LDpred-inf
m <- length(beta_gwas)
coeff <- N * h2 / m
corr2 <- corr + Matrix::Diagonal(ncol(corr), 1 / coeff)
betas_blup <- as.vector(Matrix::solve(corr2, beta_gwas))
pred <- big_prodVec(G, betas_blup, ind.row = ind.val)
plot(pred, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
cor(pred, y2[ind.val])**2  # 0.376 -> 0.38


# LDpred-fast
corr2 <- as(corr, "dgTMatrix")
corr_as_list <- split(data.frame(i = corr2@i + 1L, r = corr2@x),
                      factor(corr2@j, ordered = TRUE))

sapply(setNames(nm = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3)), function(p) {

  print(p)

  C1 <- (coeff + 1) / (coeff + p)
  C2 <- coeff / (coeff + 1)
  V_nc <- C2 / N * (1 - C2 * h2)
  V_c <- C2 * h2 / (m * p) + V_nc
  # Should use p-transformed betas instead?
  d_beta_ncaus <- dnorm(betas_blup, sd = sqrt(V_nc))
  d_beta_caus <- dnorm(betas_blup, sd = sqrt(V_c))
  d_caus_beta <- d_beta_caus * p / (d_beta_caus * p + d_beta_ncaus * (1 - p))
  beta_jms <- betas_blup * C1 * d_caus_beta

  pred5 <- big_prodVec(G, beta_jms, ind.row = ind.val)
  plot(pred5, y2[ind.val], pch = 20); abline(0, 1, col = "red", lwd = 3)
  round(100 * cor(pred5, y2[ind.val])**2, 1)
})

