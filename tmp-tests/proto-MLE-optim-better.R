library(bigsnpr)
bigsnp <- snp_attachExtdata()
G <- bigsnp$genotypes

simu <- snp_simuPheno(G, 0.2, 500, alpha = 0.5)

log_var <- log(big_colstats(G, ind.col = simu$set)$var)
beta2 <- simu$effects^2

FUN <- function(x, log_var, beta2) {
  S <- 1 + x[[1]]; sigma2 <- x[[2]]
  S * sum(log_var) + length(log_var) * log(sigma2) + sum(beta2 / exp(S * log_var)) / sigma2
}

DER <- function(x, log_var, beta2) {
  S <- 1 + x[[1]]; sigma2 <- x[[2]]
  res1 <- sum(log_var) - sum(log_var * beta2 / exp(S * log_var)) / sigma2
  res2 <- length(log_var) / sigma2 - sum(beta2 / exp(S * log_var)) / sigma2^2
  c(res1, res2)
}


optim(par = c(-0.5, 0.2 / 500), fn = FUN, method = "L-BFGS-B",
      lower = c(-2, 0.2 / 5000), upper = c(1, 0.2 / 50),
      log_var = log_var, beta2 = beta2)

optim(par = c(-0.5, 0.2 / 500), fn = FUN, method = "L-BFGS-B", gr = DER,
      lower = c(-2, 0.2 / 5000), upper = c(1, 0.2 / 50),
      log_var = log_var, beta2 = beta2)  # this one is best

optim(par = c(-0.5, 0.2 / 500), fn = FUN, gr = DER,
      # lower = c(-2, 0.2 / 5000), upper = c(1, 0.2 / 50),
      log_var = log_var, beta2 = beta2)

Rcpp::sourceCpp('tmp-tests/proto-MLE-optim-C++.cpp')
test_MLE(log_var, beta2, c(-1, 0.2 / 500),
         lower = c(-2, 0.2 / 5000), upper = c(1, 0.2 / 50))


microbenchmark::microbenchmark(
  R = optim(par = c(-1, 0.2 / 500), fn = FUN, method = "L-BFGS-B", gr = DER,
        lower = c(-2, 0.2 / 5000), upper = c(1, 0.2 / 50),
        log_var = log_var, beta2 = beta2)$par,
  C = test_MLE(log_var, beta2, c(-1, 0.2 / 500),
           lower = c(-2, 0.2 / 5000), upper = c(1, 0.2 / 50))
)

### without caching:
# Unit: microseconds
# expr    min      lq     mean  median      uq    max neval
#    R 1483.5 1580.90 1753.977 1643.35 1772.05 6920.2   100
#    C  902.8  928.05  980.327  943.30  987.50 1378.7   100
### caching does not really help because the other sums take much more time
### to compute (because of the exp)


Rcpp::sourceCpp('src/ldpred2-auto.cpp')
x <- c(0, 0.2 / 500); MLE_alpha(x, 0:499, log_var, simu$effects, FALSE)
x
