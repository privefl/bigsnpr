library(bigsnpr)

test <- AttachBigSNP("test_doc")


simu.qt <- function(X, ind, effects, h) {
  X.ind <- X[, ind]
  p.ind <- colMeans(X.ind, na.rm = T) / 2
  mean <- 2*p.ind
  sd <- sqrt(2*p.ind*(1-p.ind))

  tmp <- sweep(sweep(X.ind, 2, mean, '-'), 2, sd, '/')
  pheno.qt <- rowSums(sweep(tmp, 2, effects, '*'), na.rm = T)

  coef.h <- 1/h^2 - 1
  return(rnorm(length(pheno.qt), pheno.qt, sqrt(abs(pheno.qt) * coef.h)))
}

X <- test$genotypes
ind <- sort(sample(ncol(X), 5))
y <- sign(simu.qt(X, ind, rep(0.5,5), 0.8))
R2 <- RsqClass(X, y)
plot(R2, type = 'h')
points(ind, R2[ind], col = "red", pch = 19)

n <- nrow(X)
for (i in 1:10) {
  abline(h = qchisq(10^(-i), 1, lower.tail = FALSE) / n,
         col = i+1, lty = 2, lwd = 2)
}

