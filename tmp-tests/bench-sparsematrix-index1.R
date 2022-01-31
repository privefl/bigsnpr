library(Matrix)
m <- 5000
corr <- as(rsparsematrix(m, m, 0.1, symmetric = TRUE), "dgTMatrix")
i <- corr@i
j <- corr@j
keep <- (i <= j)
i <- i[keep]
j <- j[keep]
x <- corr@x[keep]
i2 <- i + 1L
j2 <- j + 1L

microbenchmark::microbenchmark(
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(m, m), symmetric = TRUE,
                       index1 = FALSE),
  Matrix::sparseMatrix(i = i2, j = j2, x = x, dims = c(m, m), symmetric = TRUE,
                       index1 = TRUE),
  times = 5
)
