library(bigsnpr)

chr22 <- snp_attach("../Dubois2010_data/celiac_chr22.rds")
G <- chr22$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
dim(G)
cov <- runonce::save_run(cov(G[]), file = "tmp-data/cov_chr22.rds")

i <- 3000
ind <- setdiff(which(cov[, i] ** 2 > 0.001), 1)


id_sub <- c(i, ind)
cov_sub <- as.matrix(cov[id_sub, id_sub])
rho <- 0.0001
glasso0 <- glassoFast::glassoFast(cov_sub, rho = rho)
str(glasso0)
W0 <- W0.2 <- glasso0$w; diag(W0) <- diag(W0) - rho
all.equal(W0, cov_sub) # Mean relative difference: 0.07337141
plot(W0, cov_sub); abline(0, 1, col = "red")
W0[1:5, 1:5]
cov_sub[1:5, 1:5]
all.equal(glasso0$wi %*% cov_sub, diag(ncol(cov_sub)))
inv <- solve(cov_sub + rho * diag(ncol(cov_sub)))
all.equal(glasso0$wi, inv)
glasso0$wi[1:5, 1:5]
inv[1:5, 1:5]
plot(diag(glasso0$wi), diag(inv)); abline(0, 1, col = "red")

inv2 <- solve(cov_sub + 10 * rho * diag(ncol(cov_sub)))
plot(diag(inv2), diag(inv)); abline(0, 1, col = "red")
glasso0.2 <- glassoFast::glassoFast(cov_sub, rho = 10 * rho)
plot(diag(glasso0.2$wi), diag(glasso0$wi)); abline(0, 1, col = "red")

diag(glassoFast::glassoFast(cov_sub + rho * diag(ncol(cov_sub)), rho = rho)$wi)
diag(glassoFast::glassoFast(cov_sub + 10 * rho * diag(ncol(cov_sub)), rho = 10 * rho)$wi)


all.equal(inv %*% cov_sub, diag(ncol(cov_sub)))
all.equal(inv2 %*% cov_sub, diag(ncol(cov_sub)))
glasso0 <- glassoFast::glassoFast(cov_sub, rho = rho)
inv_inv <- solve(glasso0$wi)
plot(inv_inv, cov_sub); abline(0, 1, col = "red")
hist(cov_sub - inv_inv)
hist(W0 - inv_inv)
hist(W0.2 - inv_inv)
inv_inv[1:5, 1:5]
W0.2[1:5, 1:5]
glasso_inv <- glassoFast::glassoFast(glasso0$wi, rho = 10 * rho)
all.equal(glasso_inv$wi, cov_sub)
sum(glasso_inv$wi == 0)
glasso_inv2 <- glassoFast::glassoFast(inv, rho = 10 * rho)
all.equal(glasso_inv2$wi, cov_sub)
sum(glasso_inv2$wi == 0)

Rcpp::sourceCpp("tmp-tests/test-glasso2.cpp")
glasso2 <- glasso(as.matrix(cov_sub), lambda = rho, 200, 200, tol = 1e-4, verbose = TRUE)
W <- glasso2[[1]]
all.equal(W, cov_sub)  # 0.07335971
all.equal(W, W0)       # 1e-4

X <- glasso2[[2]]

plot(X, glasso0$wi)
tmp <- 1 / (1 + rho - colSums(X * (W + rho * diag(ncol(W)))))
plot(tmp, diag(glasso0$wi))
all.equal(tmp, diag(glasso0$wi))  # 0.0002222643
X2 <- sweep(X, 2, -tmp, '*'); diag(X2) <- tmp
X3 <- (X2 + t(X2)) / 2
all.equal(X3, glasso0$wi)  # 0.001191556
plot(X3, glasso0$wi); abline(0, 1, col = "red")
all.equal(X3 %*% cov_sub, diag(ncol(cov_sub))) # 0.841867

microbenchmark::microbenchmark(
  glassoFast::glassoFast(cov_sub, rho = 0.01),
  glasso(as.matrix(cov_sub), lambda = 0.01, 200, 200, tol = 1e-4, verbose = FALSE)
)



