bdiag_sym <- function(list_corr) {
  num_tot <- sum(sapply(list_corr, function(.) length(.@x)))
  I <- integer(num_tot)
  X <- double(num_tot)
  P <- list()
  offset_I <- 0L
  offset_P <- 0
  for (k in seq_along(list_corr)) {
    corr_k <- list_corr[[k]]
    len <- length(corr_k@x)
    ind <- (1 + offset_P):(len + offset_P)
    I[ind] <- corr_k@i + offset_I
    X[ind] <- corr_k@x
    P[[k]] <- utils::head(corr_k@p, -1) + offset_P
    # print(typeof(corr_k@i + offset_I))
    offset_I <- offset_I + nrow(corr_k)
    offset_P <- offset_P + len
  }
  P[[k + 1]] <- num_tot
  Matrix::sparseMatrix(i = I, p = unlist(P), x = X, symmetric = TRUE,
                       check = FALSE, index1 = FALSE)
}

N <- sample(2e4 + -10:10, size = 10, replace = TRUE)
list_corr <- lapply(N, function(n) {
  Matrix::rsparsematrix(n, n, density = 0.002, symmetric = TRUE)
})

system.time(test <- bdiag_sym(list_corr))
system.time(true <- Matrix::bdiag(list_corr))

all.equal(test, true)
