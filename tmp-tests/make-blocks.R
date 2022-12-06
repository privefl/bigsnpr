corr0 <- readRDS("~/../Documents/LD_with_blocks_chr2.rds")

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
  ind <- which(corr[, j]^2 > 0.1)

  repeat {
    ld_with_ind <- Matrix::rowSums(corr[, ind, drop = FALSE]^2)
    thr <- 1.2 + length(ind) * r2_noise
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
