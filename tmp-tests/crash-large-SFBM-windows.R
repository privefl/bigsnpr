corr0 <- readRDS("~/../Downloads/LD_REF/ld-ref/LD_chr1.rds")

library(bigsnpr)
gc(); file.remove("tmp-data/very_large_corr.sbk")
corr <- as_SFBM(corr0, "tmp-data/very_large_corr", compact = TRUE)

for (k in seq_len(20)) {
  print(k)  # was crashing at k=3 before on Windows
  print(system.time(corr$add_columns(corr0, nrow(corr))))
}

log2(corr$nval)
file.size(corr$backingfile) / 1024^3



# NO problem with bdiag_sym() -> so a problem with appending?

bdiag_sym <- function(list_corr) {

  num_tot <- sum(sapply(list_corr, function(.) length(.@x)))
  I <- integer(num_tot)
  X <- double(num_tot)
  P <- list()
  offset_I <- 0L
  offset_P <- 0L

  for (k in seq_along(list_corr)) {
    corr_k <- list_corr[[k]]
    len <- length(corr_k@x)
    ind <- (1 + offset_P):(len + offset_P)
    I[ind] <- corr_k@i + offset_I
    X[ind] <- corr_k@x
    P[[k]] <- utils::head(corr_k@p, -1) + offset_P
    offset_I <- offset_I + nrow(corr_k)
    offset_P <- offset_P + len
  }
  P[[k + 1]] <- num_tot

  corr <- new("dsCMatrix", uplo = "U")
  dim <- sum(sapply(list_corr, ncol))
  corr@Dim <- c(dim, dim)
  corr@i <- I
  corr@p <- unlist(P)
  corr@x <- X

  corr
}

corr2 <- bdiag_sym(list(corr0, corr0, corr0, corr0, corr0))
corr <- as_SFBM(corr2, "tmp-data/very_large_corr", compact = TRUE)
