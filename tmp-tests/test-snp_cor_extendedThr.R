library(testthat)
library(bigsnpr)
bigsnp <- snp_attachExtdata()
G <- bigsnp$genotypes
POS <- bigsnp$map$physical.pos
corr0 <- big_cor(G)

thr_r2 <- THR <- 0.02

bigassertr::message2("Initial computation of all pairwise r2 > %s..", thr_r2)
corr <- snp_cor(G, thr_r2 = THR, infos.pos = POS, size = 500, ncores = 2)

ind <- Matrix::which(corr != 0, arr.ind = TRUE)
prev_list_keep <- split(ind[, "row"], factor(cols_along(corr))[ind[, "col"]])
rm(ind)

repeat {
  corr_T <- as(Matrix::drop0(corr, tol = sqrt(thr_r2)), "generalMatrix")
  new_list_keep <- bigsnpr:::find_indirect_corr(corr_T@p, corr_T@i, ncores = 2)
  # user  system elapsed
  # 2.69    0.00    1.38

  overlap <- Matrix::which((corr_T %*% corr_T) != 0, arr.ind = TRUE)
  check_list <- split(overlap[, "row"], factor(cols_along(corr))[overlap[, "col"]])
  print(c(nrow(overlap), sum(lengths(new_list_keep))))
  print(all.equal(new_list_keep, unname(check_list)))

  for (k in seq_along(prev_list_keep)) {
    diff_ind <- setdiff(new_list_keep[[k]], prev_list_keep[[k]])
    prev_list_keep[[k]] <- new_list_keep[[k]]
    new_list_keep[[k]] <- diff_ind
  }

  nb_to_add <- sum(lengths(new_list_keep))
  if (nb_to_add == 0) break

  bigassertr::message2("Computing %s additional pairwise correlations..", nb_to_add)
  corr <- corr + bigsnpr:::snp_corInd(G, list_ind = new_list_keep, ncores = 2)
}
# Computing 365306 additional pairwise correlations..
# Computing 257584 additional pairwise correlations..
# Computing 542262 additional pairwise correlations..
# Computing 798160 additional pairwise correlations..
# Computing 556546 additional pairwise correlations..
# Computing 208962 additional pairwise correlations..
# Computing 76808 additional pairwise correlations..
# Computing 32964 additional pairwise correlations..
# Computing 15656 additional pairwise correlations..
# Computing 7278 additional pairwise correlations..
# Computing 2506 additional pairwise correlations..
# Computing 994 additional pairwise correlations..
# Computing 342 additional pairwise correlations..

test <- snp_cor_extendedThr(G, thr_r2 = THR, infos.pos = POS, size = 500, ncores = 2)
all.equal(test, corr)

corr_T <- as(Matrix::drop0(test, tol = sqrt(THR)), "generalMatrix")
overlap <- Matrix::which((corr_T %*% corr_T) != 0, arr.ind = TRUE)
check_list <- split(overlap[, "row"], factor(cols_along(corr))[overlap[, "col"]])
print(all.equal(prev_list_keep, check_list))
expect_equal(test[overlap], corr0[overlap])
# expect_true(all(test[overlap] != 0))

ind <- overlap[test[overlap] == 0, ]

corr0[ind]
test[ind]
which(test[, ind[1, "row"]]**2 > THR & test[, ind[1, "col"]]**2 > THR)
c(test[4392, 1429]^2, test[4392, 48]^2, test[1429, 48]^2, test[48, 1429]^2)

bigsnpr:::any_overlap(corr_T@p, corr_T@i, 4392 - 1, 1429 - 1)
new_list_keep <- bigsnpr:::find_indirect_corr(corr_T@p, corr_T@i, ncores = 4)
new_list_keep[4392]

ind <- Matrix::which(test != 0, arr.ind = TRUE)
prev_list_keep <- split(ind[, "row"], factor(cols_along(test))[ind[, "col"]])
rm(ind)

k <- 4392
setdiff(new_list_keep[[k]], prev_list_keep[[k]])

get_list <- function(mat) {
  ind <- Matrix::which(mat != 0, arr.ind = TRUE)
  unname(split(ind[, "row"], factor(cols_along(test))[ind[, "col"]]))
}

microbenchmark::microbenchmark(
  bigsnpr:::find_indirect_corr(corr_T@p, corr_T@i, ncores = 2),
  get_list(corr_T %*% corr_T),
  get_list(Matrix::crossprod(corr_T)),
  get_list(Matrix::tcrossprod(corr_T)),
  times = 5,
  check = "equal"
)
# Unit: milliseconds
#                                                         expr       min        lq
# bigsnpr:::find_indirect_corr(corr_T@p, corr_T@i, ncores = 2) 1252.3251 1271.7342
#                                  get_list(corr_T %*% corr_T)  522.6025  530.8265
#                          get_list(Matrix::crossprod(corr_T))  433.8930  480.2534
#                         get_list(Matrix::tcrossprod(corr_T))  445.6235  463.8032
#      mean    median        uq       max neval
# 1319.1064 1321.6068 1329.5871 1420.2786     5
#  565.0707  533.8133  578.5574  659.5537     5
#  497.2513  486.2713  537.1000  548.7389     5
#  532.3320  493.2715  498.5039  760.4579     5


THR <- runif(1, 0.01, 0.5)
SIZE <- round(runif(1, 100, 2000))
ind.row <- sample(nrow(G), 600, replace = TRUE)
ind.col <- sort(sample(ncol(G), 1000, replace = TRUE))
test <- snp_cor_extendedThr(G, thr_r2 = THR, infos.pos = POS[ind.col], size = SIZE,
                            ind.row = ind.row, ind.col = ind.col, ncores = NCORES)

corr_T <- as(Matrix::drop0(test, tol = sqrt(THR)), "generalMatrix")
list_keep <- bigsnpr:::find_indirect_corr(corr_T@p, corr_T@i, ncores = NCORES)
overlap <- Matrix::which((corr_T %*% corr_T) != 0, arr.ind = TRUE)
check_list <- split(overlap[, "row"], factor(cols_along(G))[overlap[, "col"]])
expect_equal(list_keep, unname(check_list))
expect_equal(test[overlap], corr0[overlap])
NULL
