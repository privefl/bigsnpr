corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))
corr <- Matrix::tril(corr)
m <- ncol(corr)

library(magrittr)
L <- get_L(corr@p, corr@i, corr@x, thr_r2 = 0.01) %>%
  { Matrix::sparseMatrix(i = .$i, j = .$j, x = .$x, dims = c(m, m + 1),
                         triangular = FALSE, index1 = FALSE) }

path <- get_C(L, min_size = 10, max_size = 50, K = 50)

path$C[1:3, ]
dim(path$C)
plot(path$C[1, ], log = "y")
