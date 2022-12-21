files <- gtools::mixedsort(
  list.files("../datasets/EUR", full.names = TRUE))
length(files)  # 1361 -> number of blocks

snp_list <- gtools::mixedsort(
  list.files("../datasets/snplist", full.names = TRUE))

ic <- 100
prec_info <- bigreadr::fread2(files[ic])
snp_info <- bigreadr::fread2(snp_list[ic])
snp_info2 <- snp_info[!duplicated(snp_info$index), ]
prec <- Matrix::sparseMatrix(i = prec_info$V1, j = prec_info$V2, x = prec_info$V3,
                             dims = rep(nrow(snp_info2), 2), index1 = FALSE,
                             symmetric = TRUE)
prec[1:5, 1:5]

keep <- Matrix::diag(prec) != 0
af <- snp_info2$EUR[!keep]
summary(maf <- pmin(af, 1 - af))  # < 0.01
af2 <- snp_info2$EUR[keep]
summary(maf2 <- pmin(af2, 1 - af2))  # > 0.01
hist(maf2, "FD")

prec2 <- prec[keep, keep]
prec2[1:5, 1:5]
Matrix::mean(prec2  == 0)  # 99.7%
dim(prec2)


ic <- order(file.size(files), decreasing = TRUE)[1]
prec_info <- bigreadr::fread2(files[ic])
prec <- Matrix::sparseMatrix(i = prec_info$V1, j = prec_info$V2, x = prec_info$V3,
                             index1 = FALSE, symmetric = TRUE)
prec[1:5, 1:5]

keep <- Matrix::diag(prec) != 0
prec2 <- prec[keep, keep]
dim(prec2)  # 11497 x 11497

Matrix::mean(prec2 == 0)  # 99.8%
(objsize <- object.size(prec2)) / 1024^2  # 1.8 Mb
RSpectra::eigs(prec2, k = 4, sigma = -0.5)$values
# 0.003378565 0.002961802 0.002528332 0.001452597

id <- rep(0, ncol(prec2))
id[[100]] <- 1
id2 <- as.matrix(id)

prec2.2 <- as(prec2, "dgCMatrix")
prec3 <- bigsparser::as_SFBM(prec2)
file.size(prec3$backingfile) / 1024^2  # 4.5 Mb
chol_prec <- Matrix::Cholesky(prec2)

object.size(chol_prec) / objsize # x23.2

chol_prec2 <- SparseM::chol(prec2.2)
object.size(chol_prec2) / objsize # x40.2


microbenchmark::microbenchmark(
  Matrix::solve(prec2, id)[, 1],
  Matrix::solve(prec2, id2)[, 1],
  # Matrix::solve(prec2.2, id)[, 1],
  Matrix::solve(prec2.2, id2)[, 1],
  Matrix::solve(chol_prec, id2)[, 1],
  bigsparser::sp_solve_sym(prec3, id),
  times = 5,
  control = list(warmup = 1),
  check = "equal"
)
# Unit: milliseconds
#                          expr       min        lq       mean    median        uq       max neval
# Matrix::solve(prec2, id)[, 1]         35.2030   35.3345   36.46586   36.3964   36.8277   38.5677     5
# Matrix::solve(prec2, id2)[, 1]        35.1693   35.5749   35.83338   36.0704   36.0794   36.2729     5
# Matrix::solve(prec2.2, id2)[, 1]      16.7371   16.7418   16.85358   16.8162   16.9660   17.0068     5
# Matrix::solve(chol_prec, id2)[, 1]    23.2959   23.6061   23.84420   23.6426   23.7059   24.9705     5
# bigsparser::sp_solve_sym(prec3, id) 1025.0267 1036.3626 1079.50496 1046.0207 1093.6912 1196.4236     5

microbenchmark::microbenchmark(
  Matrix::solve(prec2, id)[, 1],
  Matrix::solve(prec2, id2)[, 1],
  # Matrix::solve(prec2.2, id)[, 1],
  Matrix::solve(prec2.2, id2)[, 1],
  Matrix::solve(chol_prec, id2)[, 1],
  bigsparser::sp_solve_sym(prec3, id, tol = 1e-6),
  times = 5,
  control = list(warmup = 1)
)
# Unit: milliseconds
#                          expr                         min       lq      mean   median       uq      max neval
# Matrix::solve(prec2, id)[, 1]                     34.9562  35.0441  35.55904  35.5757  35.7243  36.4949     5
# Matrix::solve(prec2, id2)[, 1]                    35.0120  35.2068  35.70062  35.7416  36.1050  36.4377     5
# Matrix::solve(prec2.2, id2)[, 1]                  16.6801  16.7488  18.30860  16.9694  18.9404  22.2043     5
# Matrix::solve(chol_prec, id2)[, 1]                23.2871  23.5515  23.93744  24.1600  24.2194  24.4692     5
# bigsparser::sp_solve_sym(prec3, id, tol = 1e-06) 585.7964 590.8296 599.02456 593.0152 612.0368 613.4448     5


object.size(prec2)           / objsize # x24.2
object.size(prec2.2)         / objsize # x145.8
object.size(chol_prec)       / objsize # x23.2
file.size(prec3$backingfile) / objsize # 2.5

Rcpp::sourceCpp("tmp-tests/test-solver.cpp")
microbenchmark::microbenchmark(
  FULL = bigsparser::sp_solve_sym(prec3, id, tol = 1e-10),
  OVERHEAD = sp_solve_sym_eigen(prec3, id, rep(0, length(id)), 1e-10, 1e4),
  times = 5,
  control = list(warmup = 1)
)
# Unit: microseconds
#     expr       min        lq      mean  median        uq       max neval
#     FULL 1030851.4 1053448.3 1062820.9 1059833 1078782.4 1091189.1     5
# OVERHEAD      72.7      86.7     416.9     179     244.4    1501.7     5

microbenchmark::microbenchmark(
  MATRIX_FREE = bigsparser::sp_solve_sym(prec3, id, tol = 1e-10),
  SP_MATRIX = sp_solve_sym_eigen2(prec2.2, id, 1e-10, 1e4),
  SP_VECTOR = sp_solve_sym_eigen3(prec2.2, id, 1e-10, 1e4),
  times = 5,
  control = list(warmup = 1),
  check = "equal"
)
# Unit: seconds
#        expr      min       lq     mean   median       uq      max neval
# MATRIX_FREE 1.021656 1.032963 1.044042 1.042033 1.043382 1.080178     5
#   SP_MATRIX 1.117750 1.128939 1.142985 1.135291 1.141687 1.191257     5
#   SP_VECTOR 1.104941 1.113508 1.129224 1.133755 1.135571 1.158347     5

microbenchmark::microbenchmark(
  MATRIX_FREE = bigsparser::sp_solve_sym(prec3, id, tol = 1e-10),
  CHOLESKY = sp_solve_sym_eigen4(prec2.2, id),
  times = 5,
  control = list(warmup = 1),
  check = "equal"
)
# Unit: seconds
#        expr      min       lq     mean   median       uq      max neval
# MATRIX_FREE 1.022654 1.038280 1.052268 1.038788 1.068432 1.093187     5
#    CHOLESKY 1.218309 1.223111 1.233765 1.226502 1.243736 1.257166     5

system.time(sp_solve_sym_eigen4(prec2.2, id, 1000))  # 7 sec
(7 - 1.23) / 999  # 0.0058

# Let's say it takes 7 sec per block
# Then it takes about 2h per iteration...
