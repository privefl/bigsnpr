A <- diag(4)
A[1, 2] <- A[2, 1] <- 0.8
A[3, 4] <- A[4, 3] <- -0.5
A2 <- A
A[2, 3] <- A[3, 2] <- 0.3
A
A2
ld <- colSums(A**2)

x0 <- rnorm(4)
b <- A %*% x0 + rnorm(1) / 1000

cbind(x0, solve(A, b), solve(A2, b))

gauss_seidel <- function(A, b, divide_diag = TRUE, niter = 20) {
  K <- length(b)
  x <- rep(0, K)
  for (iter in 1:niter) {
    # print(x)
    for (j in 1:K) {
      x[j] <- (b[j] - A[j, -j] %*% x[-j]) / `if`(divide_diag, A[j, j], 1)
    }
    # c(crossprod(x, A %*% x), crossprod(x, b))
  }
  x
}

cbind(x0, solve(A, b), solve(A2, b), gauss_seidel(A, b), gauss_seidel(A2, b))

solve(A[1:3, 1:3], b[1:3])
solve(A[2:4, 2:4], b[2:4])

K <- length(b)
x <- rep(0, K)
for (iter in 1:30) {
  x1 <- x[1:3]; x2 <- x[2:4]
  # print(x)
  for (j in 1:3) x1[j] <- (b[1:3][j] - A[1:3, 1:3][j, -j] %*% x1[-j])
  for (j in 1:3) x2[j] <- (b[2:4][j] - A[2:4, 2:4][j, -j] %*% x2[-j])
  # w2 <- c(crossprod(A[c(1, 3), 2]), crossprod(A[c(3, 4), 2]))
  w2 <- c(crossprod(A[c(1, 3), 2], solve(A[c(1, 3), c(1, 3)], A[c(1, 3), 2])),
          crossprod(A[c(4, 3), 2], solve(A[c(4, 3), c(4, 3)], A[c(4, 3), 2])))
  # w3 <- c(crossprod(A[c(1, 2), 3]), crossprod(A[c(2, 4), 3]))
  w3 <- c(crossprod(A[c(1, 2), 3], solve(A[c(1, 2), c(1, 2)], A[c(1, 2), 3])),
          crossprod(A[c(4, 2), 3], solve(A[c(4, 2), c(4, 2)], A[c(4, 2), 3])))
  w2 <- w2 / sum(w2); w3 <- w3 / sum(w3)
  x <- c(x1[1], x1[2] * w2[1] + x2[1] * w2[2], x1[3] * w3[1] + x2[2] * w3[2], x2[3])
  print(cbind(c(x1, 0), c(0, x2), x))
}
(res <- cbind(x0, A = solve(A, b), A2 = solve(A2, b),
              Al = solve(A + 0.01 * diag(4), b),
              As = solve(0.9 * A + 0.1 * diag(4), b), x))

FCR <- 45
curve(FCR + (185 - FCR) * x, from = 0.5); curve(185 * x, add = TRUE, col = "blue")
62 + (185 - 62) * c(0.6, 0.7)
