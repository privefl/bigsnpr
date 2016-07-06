# Simulating some data
N1 <- 3000
N2 <- 1000
x1 <- rnorm(N1, 0.5)
x2 <- rnorm(N2, -0.5)
x <- c(x1, x2)
x.big <- as.big.matrix(cbind(x, x+1, 2*x))

# In a case of classification
y <- c(rep(1, N1), rep(-1, N2))
res <- matrix(0, 3, ncol(x.big))
for (j in 1:ncol(x.big)) {
  mylm <- lm(y ~ x.big[, j],
             weights = c(rep(N2, N1), rep(N1, N2)))
  res[1:2, j] <- mylm$coefficients
  res[3, j] <- summary(mylm)$r.squared
}
print(res)
print(rbind(CoeffsClass(x.big, y), RsqClass(x.big, y)))

# In a case of regression
y2 <- x + rnorm(length(x), 0, 0.1)
res <- matrix(0, 3, ncol(x.big))
for (j in 1:ncol(x.big)) {
  mylm <- lm(y2 ~ x.big[, j])
  res[1:2, j] <- mylm$coefficients
  res[3, j] <- summary(mylm)$r.squared
}
print(res)
print(rbind(CoeffsReg(x.big, y2), RsqReg(x.big, y2)))

# With only half of the data
ind.train <- sort(sample(length(x), length(x) / 2))
res <- matrix(0, 3, ncol(x.big))
for (j in 1:ncol(x.big)) {
  mylm <- lm(y2[ind.train] ~ x.big[ind.train, j])
  res[1:2, j] <- mylm$coefficients
  res[3, j] <- summary(mylm)$r.squared
}
print(res)
print(rbind(CoeffsReg(x.big, y2, ind.train),
            RsqReg(x.big, y2, ind.train)))
