bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
obj.bed <- bed(bedfile)
ind.row <- rows_along(obj.bed)
ind.col <- cols_along(obj.bed)
stats <- bigsnpr:::bed_stats(obj.bed, ind.row, ind.col)
center <- stats$sum / stats$nb_nona_col
scale <- sqrt(stats$var * (stats$nb_nona_col - 1)) *
  length(ind.row) / stats$nb_nona_col

r1 <- bigsnpr:::cpMatVec4(obj.bed, ind.row, ind.col, center, scale, U[, 1])
all.equal(r1 * sqrt(nrow(G) - 2), big_univLinReg(G, U[, 1])$score)
r1_2 <- r1 / sqrt(1 - r1^2)
all.equal(r1_2 * sqrt(nrow(G) - 2), big_univLinReg(G, U[, 1])$score)

NCORES <- 2
blocks <- bigsnpr:::CutBySize(length(ind.col), nb = NCORES)
library(doParallel)
registerDoParallel(cl <- makeCluster(NCORES))
r_all <- foreach(ic = rows_along(blocks), .combine = "rbind") %:%
  foreach(jc = cols_along(U), .combine = "cbind") %dopar% {
    ind <- bigsnpr:::seq2(blocks[ic, ])
    r1 <- bigsnpr:::cpMatVec4(obj.bed, ind.row, ind.col[ind],
                              center[ind], scale[ind], U[, jc])
    r1_2 <- r1 / sqrt(1 - r1^2)
  }
stopCluster(cl)
all.equal(r_all[, 1] * sqrt(nrow(G) - 2), big_univLinReg(G, U[, 1])$score)
