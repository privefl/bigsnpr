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

require(foreach)
tmp <- foreach(i = 1:100, .combine = 'rbind') %do% {
  RsqClass(X, y, sort(sample(n, replace = TRUE)))
}

R2.min <- apply(tmp, 2, min)
plot(R2.min, type = 'h')
points(ind, R2.min[ind], col = "red", pch = 19)

R2.mean <- apply(tmp, 2, mean)
plot(R2.mean, type = 'h')
points(ind, R2.mean[ind], col = "red", pch = 19)

R2.median <- apply(tmp, 2, median)
plot(R2.median, type = 'h')
points(ind, R2.median[ind], col = "red", pch = 19)

R2.q5 <- apply(tmp, 2, function(x) quantile(x, 0.05))
plot(R2.q5, type = 'h', ylim = c(0, 0.02))
points(ind, R2.q5[ind], col = "red", pch = 19)
