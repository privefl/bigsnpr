library(bigsnpr)

load(url("https://www.dropbox.com/s/c13uygnjh6yh7vf/to-test-ldpred2.RData?raw=1"))

abs_colsums <- Matrix::colSums(abs(corr))

N <- df_beta$n_eff
scale <- sqrt(N * df_beta$beta_se^2 + df_beta$beta^2)
beta_hat <- df_beta$beta / scale

ic <- 1
new_beta <- crossprod(corr[-ic, ic], solve(corr[-ic, -ic], beta_hat[-ic]))

ind <- sample(cols_along(corr), 50)

library(future.apply)
plan("multisession", workers = 6)
new_beta <- future_sapply(ind, function(ic) {
  crossprod(corr[-ic, ic], solve(corr[-ic, -ic], beta_hat[-ic]))
})
plan("sequential")

plot(new_beta, beta_hat[ind]); abline(0, 1, col = "red", lwd = 2)
cor(new_beta, beta_hat[ind])  # 88.2%

corr2 <- corr + Matrix::Diagonal(ncol(corr), 0.1)

plan("multisession", workers = 6)
new_beta2 <- future_sapply(ind, function(ic) {
  crossprod(corr[-ic, ic], solve(corr2[-ic, -ic], beta_hat[-ic]))
})
plan("sequential")
plot(new_beta2, beta_hat[ind]); abline(0, 1, col = "red", lwd = 2)
cor(new_beta2, beta_hat[ind])  # 97.4%

corr3 <- corr + Matrix::Diagonal(ncol(corr), 0.5)

plan("multisession", workers = 6)
new_beta3 <- future_sapply(ind, function(ic) {
  crossprod(corr[-ic, ic], solve(corr3[-ic, -ic], beta_hat[-ic]))
})
plan("sequential")
plot(new_beta3, beta_hat[ind]); abline(0, 1, col = "red", lwd = 2)
cor(new_beta3, beta_hat[ind])  # 97.4%
