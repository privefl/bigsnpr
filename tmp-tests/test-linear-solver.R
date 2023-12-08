corr <- runonce::save_run(
  snp_cor(chr22$genotypes, infos.pos = POS2, size = 3 / 1000, ncores = 6),
  file = "tmp-data/corr_chr22.rds"
)

ind <- setdiff(order(corr[, 1] ** 2, decreasing = TRUE), 1)

id_sub <- 1
ic <- 0

j <- ind[ic <- ic + 1]

id_sub <- c(id_sub, j)
corr_sub <- as.matrix(corr[id_sub, id_sub])
glasso <- glassoFast::glassoFast(corr_sub, rho = 0.01)
glasso$wi

ind2 <- which(glasso$wi[, 1] != 0)
id_sub2 <- id_sub[ind2]
length(id_sub2) / length(id_sub)
corr_sub2 <- as.matrix(corr[id_sub2, id_sub2])
glasso2 <- glassoFast::glassoFast(corr_sub2, rho = 0.01)
plot(glasso2$wi[, 1], glasso$wi[ind2, ind2][, 1])
all.equal(glasso2$wi[, 1], glasso$wi[ind2, ind2][, 1])

glasso3 <- bigutilsr::regul_glasso(corr_sub2, lambda = 0.01)
all.equal(glasso3,   corr_sub2)
all.equal(glasso2$w, corr_sub2)

glasso4 <- glasso2$w; diag(glasso4) <- 1
all.equal(glasso3, glasso2$w)
all.equal(glasso3, glasso4)
glasso3[1:5, 1:5]
glasso2$w[1:5, 1:5]

microbenchmark::microbenchmark(
  glasso2 <- glassoFast::glassoFast(corr_sub2, rho = 0.01),
  glasso3 <- bigutilsr::regul_glasso(corr_sub2, lambda = 0.01)
)

corr_sub2_sfbm <- bigsnpr::as_SFBM(corr[id_sub2, id_sub2], compact = TRUE)
b <- rep(0, length(id_sub2)); b[1] <- 1
test_solve <- bigsparser::sp_solve_sym(corr_sub2_sfbm, b, add_to_diag = 0)
plot(glasso2$wi[, 1], test_solve); abline(0, 1, col = "red")

Rcpp::sourceCpp("tmp-tests/test-linear-solver2.cpp",
                showOutput = FALSE, echo = FALSE)
test_solve3 <- test_solver(as(corr[id_sub2, id_sub2], "generalMatrix"), b)

plot(glasso2$wi[, 1], test_solve3); abline(0, 1, col = "red")
plot(test_solve3, test_solve)
all.equal(test_solve3, test_solve)

corr_sub3 <- as(corr[id_sub2, id_sub2], "generalMatrix")

microbenchmark::microbenchmark(
  glasso2 <- glassoFast::glassoFast(corr_sub2, rho = 0.01),
  glasso3 <- bigutilsr::regul_glasso(corr_sub2, lambda = 0.01),
  test_solve <- bigsparser::sp_solve_sym(corr_sub2_sfbm, b),
  test_solve2 <- bigsparser::sp_solve_sym(corr_sub2_sfbm, b, add_to_diag = 0.01),
  # as(corr[id_sub2, id_sub2], "generalMatrix"),
  test_solve3 <- test_solver(corr_sub3, b)
)


corr2 <- as(corr, "generalMatrix")
b2 <- rep(0, ncol(corr2)); b2[1] <- 1
system.time(
  test_solve4 <- test_solver(corr2, b2)
)
plot(test_solve4[ind2], test_solve3)

system.time(
  test_solve5 <- test_solver(corr2 + Matrix::Diagonal(n = ncol(corr2), x = 0.01), b2)
) # 39 sec
plot(test_solve4, test_solve5)

system.time(
  test_solve6 <- test_solver(corr2 + Matrix::Diagonal(n = ncol(corr2), x = 0.1), b2)
) # 8 sec
plot(test_solve6, test_solve5); abline(0, 1, col = "red")

Rcpp::sourceCpp("tmp-tests/test-linear-solver3.cpp")
corr3 <- as.matrix(corr2)
system.time(
  test_solve7 <- test_solver_dense(corr3, b2)
) # 87 sec
all.equal(test_solve7, test_solve4)

corr4 <- corr3; diag(corr4) <- diag(corr3) + 0.01
system.time(
  test_solve8 <- test_solver_dense(corr4, b2)
) # 21 sec
all.equal(test_solve8, test_solve5)

system.time(
  test_solve9 <- test_solver_dense_sparse(corr4)
) # same time
all.equal(test_solve8, test_solve9)

system.time(
  test_solve10 <- test_ConjugateGradient(corr4)
) # 28 sec
all.equal(test_solve10, test_solve9)

system.time(
  test_solve11 <- test_BiCGSTAB(corr4)
) # 225 sec
all.equal(test_solve11, test_solve9)  # 1e-4 relative diff

system.time(
  test_solve12 <- test_GMRES(corr4)
) # 119 sec
all.equal(test_solve12, test_solve9)  # 1e-2 relative diff
