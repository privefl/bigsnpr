corr0 <- readRDS("~/../Documents/LD_with_blocks_chr22.rds")

get_mode <- function(x) {
  d <- density(log10(x))
  10^d$x[which.max(d$y)]
}

(r2_noise <- get_mode(corr0@x^2))

corr <- as(corr0, "dgCMatrix")
nval <- length(corr@x)
diag(corr) <- ncol(corr)
keep <- seq_len(ncol(corr))
ic <- 1
blocks <- list()

while (length(keep) > 0) {

  ld <- Matrix::colSums(corr^2)
  j <- which.max(ld)
  ind <- which(corr[, j]^2 > 0.2)

  repeat {
    ld_with_ind <- Matrix::rowSums(corr[, ind, drop = FALSE]^2)
    thr <- 1.5 + length(ind) * r2_noise
    new_ind <- which(ld_with_ind > thr)
    if (length(new_ind) > length(ind)) {
      ind <- new_ind
      print(length(ind))
    } else {
      blocks[[ic]] <- keep[ind]
      print(ic <- ic + 1)
      keep <- setdiff(keep, keep[ind])
      corr <- corr[-ind, -ind, drop = FALSE]
      break
    }
  }
}

length(blocks)
plot(lengths(blocks), pch = 20, log = "y")
crossprod(lengths(blocks)) / nval
# 53% for chr22
# 63% for chr9
# 95.5% for chr22


blocks2 <- blocks[order(purrr::map_int(blocks, 1))]
plot(unlist(blocks2), col = rep(seq_along(blocks2), lengths(blocks2)))

best_grp <- rep(NA, ncol(corr0))
for (ic in seq_along(blocks)) best_grp[blocks[[ic]]] <- ic
corr0T <- as(corr0, "dgTMatrix")
corr0T@x <- ifelse(best_grp[corr0T@i + 1L] == best_grp[corr0T@j + 1L], corr0T@x, 0)
corr2 <- as(Matrix::drop0(corr0T), "symmetricMatrix")
object.size(corr0) / object.size(corr2)

RSpectra::eigs(corr0, k = 2)$values
# 106.13132  83.42977
RSpectra::eigs(corr0, k = 2, sigma = -0.05)$values
# -0.007954082 -0.009787068

RSpectra::eigs(corr2, k = 2)$values
# 106.13132  83.42977
RSpectra::eigs(corr2, k = 2, sigma = -0.01)$values
# -0.001336321 -0.002949230 -> -1.168228e-17 -1.208885e-16
