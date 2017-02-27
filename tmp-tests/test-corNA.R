require(bigsnpr)

X <- big.matrix(1000, 200, type = "char")
X[] <- sample(c(0:2, NA), size = length(X), replace = TRUE)
n <- nrow(X)
m <- ncol(X)

result2 <- Hmisc::rcorr(X[,])
round(result2$r[, 50] * 100, 1)
round(result2$P[, 50] * 100, 1)
plot(result2$r[, 50], result2$P[, 50])

# sum(round(result2$P[, 50] * 100, 1) < 5, na.rm = TRUE)

test <- pchisq(result2$n * result2$r^2, df = 1, lower.tail = FALSE)
plot(test, result2$P)

q.alpha <- qchisq(1, df = 1, lower.tail = FALSE)
thrs <- q.alpha / 1:1000

Rcpp::sourceCpp('src/corr.cpp')
test <- corMat(X@address, rowInd = 1:1000, colInd = 1:200,
               size = 200, thr = thrs)
test[1:5, 1:5]

true <- cor(X[,], use = "pairwise.complete.obs")
true[1:5, 1:5]
all.equal(test[upper.tri(test)], true[upper.tri(true)])
