library(bigsnpr)
bigsnp <- snp_attachExtdata()
G <- bigsnp$genotypes

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

alpha <- 0
simu <- snp_simuPheno(G, 0.2, 500, alpha = alpha)
log_var <- log(big_colstats(G, ind.col = simu$set)$var)
beta2 <- simu$effects^2
all_est <- replicate(200, {
  ind <- sample(500, replace = TRUE)
  optim(par = c(-0.5, 0.2 / 500), fn = FUN, gr = DER, method = "L-BFGS-B",
        lower = c(-1.5, 0.2 / 5000), upper = c(0.5, 0.2 / 50),
        log_var = log_var[ind], beta2 = beta2[ind])$par
})

Rcpp::sourceCpp('src/ldpred2-auto.cpp')
all_est2 <- replicate(200, {
  MLE_alpha(par = c(-0.5, mean(all_est[2, ])), ind_causal = seq_along(log_var) - 1L,
            log_var = log_var, curr_beta = simu$effects,
            boot = TRUE)[1:2] - 1:0
})
hist(all_est[1, ], col = scales::alpha("orange", 0.3))
hist(all_est2[1, ], col = scales::alpha("blue", 0.3), add = TRUE)
hist(all_est[2, ], col = scales::alpha("orange", 0.3))
hist(all_est2[2, ], col = scales::alpha("blue", 0.3), add = TRUE)
