require(HiLMM)

data_sim <- data_simu(n = 1e3, N = 1e4, eta_star = 0.99, q = 0.01)

# print(system.time(
#   test <- estim_herit(Y = data_sim$Y, W = data_sim$W)
# ))
# str(test)


estim_herit2 <- function(Y, W, eta_init = c(0.1, 0.5, 0.9), nb_iter = 20) {
  n <- nrow(W)
  N <- ncol(W)

  Z <- scale(W, center = TRUE, scale = TRUE)
  M <- tcrossprod(Z) / N
  eigs <- eigen(M, symmetric = TRUE)
  lambda <- eigs$values
  Y_tilde <- crossprod(eigs$vectors, Y)

  tmp1 <- lambda - 1
  Y_tilde2 <- Y_tilde^2

  NR_step <- function(eta) {
    tmp2 <- (eta * tmp1 + 1)
    tmp3 <- tmp1 / tmp2
    tmp4 <- Y_tilde2 / tmp2

    A <- sum(tmp4 * tmp3) / sum(tmp4) - sum(tmp3) / n
    B <- -2 * sum(tmp4 * tmp3^2) / sum(tmp4) +
      (sum(tmp4 * tmp3) / sum(tmp4))^2 +
      sum(tmp3^2) / n

    eta - A / B
  }

  l <- length(eta_init)
  eta_vec <- numeric(l)

  for (i in 1:l) {
    eta <- eta_init[i]
    for (nb in 1:nb_iter) {
      eta <- NR_step(eta)
      #print(eta)
    }
    eta_vec[i] <- eta
  }

  ind_init <- which.min(abs(eta_vec - 0.5))
  eta_chap <- eta_vec[ind_init]
  w <- tmp1 / (eta_chap * tmp1 + 1)
  s_eta <- sqrt(2/(sum(w^2) - 1/n * sum(w)^2))
  dif <- qnorm(0.975) * s_eta
  list(heritability = eta_chap,
       CI_low = eta_chap - dif,
       CI_up = eta_chap + dif,
       standard_deviation = s_eta)
}

print(system.time(
  test2 <- estim_herit2(Y = data_sim$Y, W = data_sim$W,
                        #eta_init = 0.00001,
                        nb_iter = 20)
))
str(test2)
#print(all.equal(test, test2))
