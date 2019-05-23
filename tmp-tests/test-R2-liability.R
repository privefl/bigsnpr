library(bigsnpr)
G <- snp_attachExtdata()$genotypes

test <- replicate(10e3, {

  h2 <- 0.5; K <- 0.1; M <- ncol(G)
  set <- sort(sample(ncol(G), size = M))
  effects <- rnorm(M, sd = sqrt(h2 / M))
  y <- drop(scale(G[, set]) %*% effects)       ## G
  y2 <- y + rnorm(nrow(G), sd = sqrt(1 - h2))  ## G + E
  y3 <- as.integer(y2 > qnorm(1 - K))          ## LTM
  var(y) / var(y2)                             ## H2

  lmv <- lm(y3 ~ y)
  P <- mean(y3)
  (R20 <- var(lmv$fitted.values) / (P * (1 - P)))   ## 0.5055679
  # summary(lmv)$r.squared
  R20 <- summary(lmv)$adj.r.squared               ## 0.499598
  thd <- qnorm(1 - K)
  zv <- dnorm(thd)
  mv <- zv / K
  mvPK <- mv * (P - K) / (1 - K)
  theta <- mvPK * (mvPK - thd)
  cv <- (K * (1 - K) / zv)^2 / (P * (1 - P))
  (R2 <- R20 * cv / (1 + R20 * cv * theta))
  c(R2, AUC(y, y3))

})

plot(t(test), pch = 20)
mean(test[1, ])
