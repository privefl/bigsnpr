library(bigsnpr)
bigsnp <- snp_attachExtdata()
G <- bigsnp$genotypes
simu <- snp_simuPheno(G, h2 = 0.2, M = 100)
y <- simu$pheno * runif(1) + runif(1)
cov <- matrix(rnorm(3 * length(y)), ncol = 3)

all_coef <- t(apply(G[], 2, function(g) {
  summary(lm(y ~ g + cov))$coef["g", ]
}))

all_coef2 <- big_univLinReg(G, y, covar.train = cov)
all.equal(as.matrix(all_coef2), all_coef[, 1:3], check.attributes = FALSE)

all_beta <- apply(G[], 2, function(g) {
  X <- cbind(1, g, cov)
  C <- solve(crossprod(X))
  c((C %*% crossprod(X, y))[[2]], C[2, 2])
})
all.equal(all_coef[, "Estimate"], all_beta[1, ])
cor(all_coef[, "Std. Error"], sqrt(all_beta[2, ]))

w <- runif(length(y), min = 1, max = 5)
all_coef_w <- t(apply(G[], 2, function(g) {
  summary(lm(y ~ g + cov, weights = w))$coef["g", ]
}))
all_coef_w_scale <- t(apply(G[], 2, function(g) {
  summary(lm(y ~ g + cov, weights = 2 * w))$coef["g", ]
}))
all.equal(all_coef_w_scale, all_coef_w)

w2 <- w / max(w)
all_beta_w <- apply(G[], 2, function(g) {
  X <- cbind(1, g, cov)
  X2 <- sweep(X, 1, w2, '*')
  C <- solve(crossprod(X2, X))
  c((C %*% crossprod(X2, y))[[2]], C[2, 2])
})
all.equal(all_coef_w[, "Estimate"], all_beta_w[1, ])
cor(all_coef_w[, "Std. Error"], sqrt(all_beta_w[2, ]))
plot(all_coef_w[, "Std. Error"], sqrt(all_beta_w[2, ])); abline(0, 1, col = "red")

Rcpp::sourceCpp("tmp-tests/test-weighted-linreg.cpp")
all_coef_w2 <- big_univLinReg_w(G, y, w, covar.train = cov)
all_coef_w2$pval <- predict(all_coef_w2, log10 = FALSE)
all.equal(as.matrix(all_coef_w2), all_coef_w, check.attributes = FALSE)
