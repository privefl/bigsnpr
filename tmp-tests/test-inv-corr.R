G <- bigstatsr::big_attachExtdata()
S <- crossprod(scale(G[], center = TRUE, scale = FALSE)) / (nrow(G) - 1)
system.time(
  sol <- glassoFast::glassoFast(S, rho = 1 / sqrt(nrow(G)), trace = TRUE)
) # 10 iter in 83 sec
hist(sol$w - S)
mean(sol$wi == 0) # 98.3%

cov2cor_scaler <- tcrossprod(sqrt(diag(solve(sol$wi))))  # also 1-2 min to run
inv_corr <- as(sol$wi * cov2cor_scaler, "dsCMatrix")
object.size(sol$wi) / object.size(inv_corr) # x76.2

id <- matrix(0, ncol(inv_corr), 1)
id[[100]] <- 1

str(inv_corr)
Matrix::solve(inv_corr, id)  # a bit slow the first time
Matrix::solve(inv_corr, id)  # fast the second time
str(inv_corr)  # because now store some SPdCholesky decomp


object.size(sol$wi) / object.size(inv_corr) # x1.9 -> but also not very small anymore..

# can be made a but faster when not using a sparse matrix
inv_corr2 <- as(inv_corr, "dgCMatrix")  # not stored as symmetric

str(id2 <- as(id, "dgCMatrix"))

chol <- as(Matrix::chol(inv_corr), "dgCMatrix")
Matrix::nnzero(chol) / length(chol)
object.size(sol$wi) / object.size(chol) # x1.9


microbenchmark::microbenchmark(
  Matrix::solve(inv_corr, id)[, 1],
  Matrix::solve(inv_corr2, id)[, 1],
  Matrix::solve(inv_corr, id2)[, 1],
  Matrix::solve(inv_corr2, id2)[, 1],
  times = 10,
  check = "equal"
)
# Unit: milliseconds
#                         expr     min      lq     mean   median      uq     max
#  Matrix::solve(inv_corr, id) 60.9341 61.5256 62.50316 61.75390 62.0496 68.6997
# Matrix::solve(inv_corr2, id) 17.3243 17.4847 18.52723 18.27585 19.2379 21.0002
