library(Hmisc)

x <- matrix(NA_real_, 1000, 10)
x[] <- runif(length(x))
x[x > 0.5] <- NA
result <- rcorr(x)
# result$r[result$n<5] <- 0 # ignore less than five observations
result$r

result$P



require(bigsnpr)

popres.chr1 <- snp_attach("backingfiles/popres_sub2.bk")
X2 <- popres.chr1$genotypes
n <- nrow(X2)
m <- ncol(X2)

X2.beg <- X2[, 1:200]
result2 <- rcorr(X2.beg)
round(result2$r[, 50] * 100, 1)
round(result2$P[, 50] * 100, 1)
plot(result2$r[, 50], result2$P[, 50])

sum(round(result2$P[, 50] * 100, 1) < 5, na.rm = TRUE)

test <- pchisq(result2$n * result2$r^2, df = 1, lower.tail = FALSE)
plot(test, result2$P)

q.alpha <- qchisq(0.1, df = 1, lower.tail = FALSE)
thrs <- q.alpha / 1:1385

test <- corMat(X2@address, rowInd = 1:1385, colInd = 1:200, size = 200, thr = thrs)
test[1:5, 1:5]

true <- cor(X2.beg, use = "pairwise.complete.obs")
true[1:5, 1:5]
true[1:5, 1:5] / test[1:5, 1:5]

result2$r[1:5, 1:5]

true2 <- result2$r[1, 2]
x1 <- X2[, 1]
nbNA1 <- sum(is.na(x1))
x2 <- X2[, 2]
nbNA2 <- sum(is.na(x2))
x1.sum <- sum(x1, na.rm = TRUE)
xx1.sum <- sum(x1^2, na.rm = TRUE)
x2.sum <- sum(x2, na.rm = TRUE)
xx2.sum <- sum(x2^2, na.rm = TRUE)
xy <- x1 * x2
nbNA3 <- sum(is.na(xy))
xy.sum <- sum(xy, na.rm = TRUE)

num <- xy.sum - x1.sum * x2.sum / (n - nbNA3)
deno1 <- (xx1.sum - x1.sum * x1.sum / (n - nbNA1)) / (n - nbNA1 - 1)
print(deno1 - var(x1, na.rm = TRUE))
deno2 <- (xx2.sum - x2.sum * x2.sum / (n - nbNA2)) / (n - nbNA2 - 1)
print(deno2 - var(x2, na.rm = TRUE))
cov(x1, x2, use = "pairwise.complete.obs")

tmp <- (x1 - x1.sum / nbNA1) * (x2 - x2.sum / nbNA2)
sum(tmp, na.rm = TRUE)
