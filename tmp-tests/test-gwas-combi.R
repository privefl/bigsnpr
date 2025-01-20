library(bigsnpr)
obj <- snp_attachExtdata()
G <- obj$genotypes
set.seed(1)
y <- rnorm(nrow(G))
w <- runif(5); w <- w / sum(w)
G_w <- G[, 1:5] %*% w

sapply(1:5, function(j) {
  summary(lm(y ~ G[, j]))$coef[2, 1:3]
})

gwas <- big_univLinReg(G, y, ind.col = 1:5)
gwas$score

(Z_w <- crossprod(gwas$score, w))
summary(mylm <- lm(y ~ G_w))$coef  # 0.570 instead of 0.266
# crossprod(scale(G_w), scale(y)) / (sd(mylm$residuals) * sqrt(nrow(G)))

var <- big_colstats(G, ind.col = 1:5)$var
(beta_w <- crossprod(gwas$estim * var, w) / var(G_w))  # both 0.0739
