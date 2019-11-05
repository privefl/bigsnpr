
popsmooth <- function(pop, t1, z, v, gamma, eimax) {

  lossf <- function(pop, eimax, z, v, gamma) {
    gamma2 <- gamma / length(pop)
    tm1 <- eimax / pop[pop > 0]  # 1 / t
    term1 <- sapply(v, function(v_i) sum(1 / (tm1 + v_i)))
    e <- 1 / v + z - gamma2 * term1
    max(abs(Re(e)), abs(Im(e)))
  }

  seq_bw <- seq_log(length(pop) / 2, sqrt(length(pop)), 30)
  library(doParallel)
  registerDoParallel(cl <- makeCluster(ncores))
  on.exit(stopCluster(cl), add = TRUE)
  pop_loss <- foreach(bw = seq_bw, .combine = 'rbind') %dopar% {
    pop <- stats::ksmooth(seq_along(pop), pop, bandwidth = bw)$y
    loss <- lossf(pop, eimax, z, v, gamma)
    list(pop, loss)
  }

  pop_loss[[which.min(pop_loss[, 2]), 1]]
}

karoui.nonsp <- function(samp.eval, m, p, n) {

  gamma <- p / n
  eimax <- samp.eval[m + 1]
  eimin <- `if`(n > p, samp.eval[n], 0)

  samp.scaled <- samp.eval[(m + 1):n] / eimax
  rev <- seq(0, 1, 0.01)
  imv <- seq(0.001, 0.01, length.out = 5)
  v <- matrix(0, length(rev), length(imv))
  z <- matrix(mean(samp.scaled) + (0 + 1i), length(rev), length(imv))
  for (i in 1:length(rev)) {
    for (j in 1:length(imv)) {
      v[i, j] <- complex(real = rev[i], imaginary = imv[j])
      if (i == 1) {
        z[i, j] <- complex(real = mean(samp.scaled),
                           imaginary = 1 / Im(v[i, j]))
      }
      if (i == 2 && j > 1)
        z[i, j] <- z[i, (j - 1)]
      if (i > 2)
        z[i, j] <- z[(i - 1), j]

      flag.err <- 1
      try({
        z[i, j] <- hdpca:::invsteiltjes(v[i, j], z[i, j], samp.scaled, m, n)
        flag.err <- 0
      }, silent = TRUE)
      if (flag.err == 1) {
        try({
          z[i, j] <- complex(real = mean(samp.scaled), imaginary = 1 / Im(v[i, j]))
          z[i, j] <- hdpca:::invsteiltjes(v[i, j], z[i, j], samp.scaled, m, n)
          flag.err <- 0
        }, silent = TRUE)
      }
      if (Im(z[i, j]) == 0 || flag.err == 1)
        stop("The inverse stieltjes transform does not converge in Karoui's algorithm")
    }
  }
  z <- as.vector(z)
  v <- sapply(z, function(z_i) {
    sum(1 / (samp.scaled - z_i)) / (n - m)
  })
  t <- seq(eimin / eimax, 1, 5 * 10 ^ -3)
  ai <- seq(eimin / eimax, 1 - 2 ^ (-7), 2 ^ (-7))
  bi <- seq(eimin / eimax + 2 ^ (-7), 1, 2 ^ (-7))
  a <- c(rep(0, (length(t) + length(ai))), 1)
  A1 <- NULL
  A2 <- NULL
  A3 <- c(rep(1, (length(t) + length(ai))), 0)
  b3 <- 1
  b1 <- NULL
  b2 <- NULL
  for (i in seq_along(z)) {

    z1 <- z[i]
    v2 <- v[i]

    coef <- p / sum(n)
    coefs1 <- t / (1 + t * v2)
    coefs2 <- 1 / v2 - log((1 + v2 * bi) / (1 + v2 * ai)) / ((bi - ai) * v2^2)
    a11 <-  Re(coef * coefs1)
    a12 <-  Re(coef * coefs2)
    a21 <- -Im(coef * coefs1)
    a22 <- -Im(coef * coefs2)
    A1 <- rbind(A1, c(a11, a12, -1))
    A1 <- rbind(A1, c(a21, a22, -1))
    A2 <- rbind(A2, c(a11, a12,  1))
    A2 <- rbind(A2, c(a21, a22,  1))

    term1 <-  Re(v2) / (Mod(v2) ^ 2) + Re(z1)
    term2 <- -Im(v2) / (Mod(v2) ^ 2) + Im(z1)
    b1 <- c(b1, term1, -term2)
    b2 <- c(b2, term1, -term2)

    if (term1 < 0 || term2 > 0) break
  }
  sign <- c(rep("<=", nrow(A1)), rep(">=", nrow(A2)), rep("=", 1))

  ERROR_MSG <- "The solution for the linear programming problem is not available in Karoui's algorithm"

  flag.err <- 0
  try({
    simp <- lpSolve::lp(direction = "min", a, rbind(A1, A2, A3), sign, c(b1, b2, b3))
    flag.err <- 1
  })
  if (flag.err == 0) stop(ERROR_MSG)

  if (simp$status != 0) {
    flag.err <- 0
    try({
      simp <- boot::simplex(a, A1, b1, A2, b2, A3, b3)
      flag.err <- 2
    })
    if (flag.err == 0) stop(ERROR_MSG)
    flag.err <- ifelse(simp$solved == 1, 2, 0)
  }
  if (flag.err == 0) stop(ERROR_MSG)

  oldw <- options(warn = -1)
  t <- seq(eimin / eimax, 1, 5 * 10 ^ -3)
  t1 <- sort(c(t, ai))
  if (flag.err == 1) {
    w1 <- simp$solution[1:length(t)]
    w2 <- simp$solution[(length(t) + 1):(length(t) + length(ai))]
  } else {
    w1 <- simp$soln[1:length(t)]
    w2 <- simp$soln[(length(t) + 1):(length(t) + length(ai))]
  }
  cuml <- sapply(t1, function(t1_j) {
    a <- sum(w1[which(t <= t1_j)])
    b <- sum(w2[which(bi <= t1_j)])
    if (max(which(bi <= t1_j)) == length(bi) || max(which(bi <= t1_j)) == -Inf) {
      d <- 0
    } else {
      c <- max(which(bi <= t1_j)) + 1
      d <- w2[c] * (t1_j - ai[c]) / (bi[c] - ai[c])
    }
    a + b + d
  })
  options(oldw)

  est <- sapply((p - m):1, function(i) {
    hdpca:::quant((i - 1) / (p - m), t1, cuml)
  })
  pop <- popsmooth(est * eimax, t1, z, v, gamma, eimax)
  c(rep(pop[1], m), pop)
}

select.nspike <- function(samp.eval, p, n, n.spikes.max) {

  if (length(samp.eval) == (n - 1)) samp.eval <- c(samp.eval, 0)
  if (length(samp.eval) != n) stop("'samp.eval' must have length n or (n-1).")

  m <- n.spikes.max
  repeat {
    pop.nonsp <- karoui.nonsp(samp.eval, m, p, n)
    eval.l <- hdpca:::l.eval(samp.eval, pop.nonsp, m, p, n)
    if (min(eval.l) > 0) {
      break
    } else {
      m <- min(which(eval.l == min(eval.l))) - 1
      if (m == 0) stop("No distant spike was found.")
    }
  }

  m
}
