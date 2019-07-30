read_bed_scaled(obj.bed, 1:10, 1:10, center[1:10], scale[1:10])
ind.row <- rows_along(obj.bed)

bed_tcrossprodSelf <- function(
  obj.bed,
  ind.row = rows_along(obj.bed),
  ind.col = cols_along(obj.bed),
  block.size = block_size(length(ind.row))
) {

  stats <- bigsnpr:::bed_stats(obj.bed, ind.row, ind.col)
  af <- stats$sum / (2 * stats$nb_nona_col)
  center <- 2 * af
  scale <- sqrt(2 * af * (1 - af))

  n <- length(ind.row)
  m <- length(ind.col)

  scale2 <- scale * sqrt(stats$nb_nona_col / n)
  scale3 <- sqrt(stats$nb_nona_row / m)

  K <- FBM(n, n, init = 0)

  intervals <- bigsnpr:::CutBySize(m, block.size)

  for (j in rows_along(intervals)) {
    ind <- bigsnpr:::seq2(intervals[j, ])
    tmp <- read_bed_scaled(obj.bed, ind.row, ind.col[ind],
                           center[ind], scale2[ind])
    bigstatsr:::incrSup2(K, tcrossprod(tmp))
  }

  # Complete the lower part of the symmetric matrix
  bigstatsr:::complete2(K)
  bigstatsr:::scaleK(K, rep(0, n), rep(0, n), scale3, 0)
  structure(K, stats = stats)
}

val <- eigen(K[], symmetric = TRUE, only.values = TRUE)$values
val
rbind(obj.svd$d, sqrt(head(val, length(obj.svd$d))))
G <- K[] / length(ind.col)
ind <- which(abs(G) > 0.35, arr.ind = TRUE)
hist(G[ind[ind[, 1] < ind[, 2], ]])

ind <- which(abs(G) > 0.5, arr.ind = TRUE)
ind[ind[, 1] < ind[, 2], ]
