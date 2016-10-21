


require(bigsnpr)

x <- AttachBigSNP("test_doc")

# scaling
p <- colmeans(x$genotypes) / 2
sd <- sqrt(2 * p * (1 - p))
Y <- sweep(sweep(x$genotypes[,], 2, 2 * p, '-'),
           2, sd, '/')



# parameters
K <- 10
L <- 2 * K
m <- ncol(Y)
n <- nrow(Y)
I <- 10

# algo
G <- matrix(rnorm(n * L), n, L) # G0
lH <- list()

for (i in 1:I) {
  lH[[i]] <- crossprod(Y, G)
  if (i < I) G <- (Y %*% lH[[i]]) / m
}

H <- do.call(cbind, lH)
H.svd <- svd(H)

T.t <- Y %*% H.svd$u
T.svd <- svd(T.t)
approx <- sweep(T.svd$u[, 1:10], 2, (T.svd$d)[1:10], '*')

pca <- prcomp(Y)
true <- pca$x[, 1:10]

plot(as.numeric(true), as.numeric(approx), pch = 19)
