set.seed(1)
A <- diag(4)
A[1, 2] <- A[2, 1] <- 0.8
A[3, 4] <- A[4, 3] <- -0.5
A_block <- A
A[2, 3] <- A[3, 2] <- 0.3
A
A_block

x0 <- rnorm(4)
b <- A %*% x0 + rnorm(1) / 1000

# iterative algorithm to solve Ax=b
Gauss_Seidel <- function(A, b, niter = 20, verbose = FALSE) {
  K <- length(b)
  x <- rep(0, K)
  for (iter in 1:niter) {
    for (j in 1:K) x[j] <- (b[j] - A[j, -j] %*% x[-j]) / A[j, j]
    if (verbose) print(x)
  }
  x
}

cbind.data.frame(x0, solve(A, b), Gauss_Seidel(A, b),
                 solve(A_block, b), Gauss_Seidel(A_block, b))



# Try overlapping block strategy
b1 <- b[1:3]; b2 <- b[2:4]; A1 <- A[1:3, 1:3]; A2 <- A[2:4, 2:4]

# w2 <- c(crossprod(A1[-2, 2]), crossprod(A2[-1, 1]))
w2 <- c(crossprod(A1[-2, 2], solve(A1[-2, -2], A1[-2, 2])),
        crossprod(A2[-1, 1], solve(A2[-1, -1], A2[-1, 1])))
# w3 <- c(crossprod(A1[-3, 3]), crossprod(A2[-2, 2]))
w3 <- c(crossprod(A1[-3, 3], solve(A1[-3, -3], A1[-3, 3])),
        crossprod(A2[-2, 2], solve(A2[-2, -2], A2[-2, 2])))
w2 <- w2 / sum(w2); w3 <- w3 / sum(w3)
list(w2, w3)  # weights to merge overlapping results for several blocks

K <- length(b)
x <- rep(0, K)
for (iter in 1:30) {
  x1 <- x[1:3]; x2 <- x[2:4]
  for (j in 1:3) x1[j] <- (b1[j] - A1[j, -j] %*% x1[-j]) / A1[j, j]
  for (j in 1:3) x2[j] <- (b2[j] - A2[j, -j] %*% x2[-j]) / A2[j, j]
  x <- c(x1[1], x1[2] * w2[1] + x2[1] * w2[2], x1[3] * w3[1] + x2[2] * w3[2], x2[3])
  print(cbind.data.frame(x1 = c(x1, 0), x2 = c(0, x2), x))
}
cbind.data.frame(x0, solve(A, b), solve(A_block, b), x)


cbind(A, bigutilsr::regul_glasso(A, 0.01))
bigutilsr::regul_glasso(A_block, 0.01)
A_block
