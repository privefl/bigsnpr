linearShrinkLWEst <- function(dat) {

  # get the number of variables and observations
  p <- ncol(dat)
  n <- nrow(dat)
  dat <- as.matrix(dat)

  # compute the sample covariance matrix
  cov_mat <- cov(dat)

  # estimate the scalers
  m_n <- mean(diag(cov_mat))
  cov2.norm2 <- norm(cov_mat, type = "F")^2
  d_n_2 <- cov2.norm2 / p - m_n^2

  dat_centered <- scale(dat, center = TRUE, scale = FALSE)
  b_n_2 <- (sum(crossprod(dat_centered^2)) - n * cov2.norm2) / (n^2 * p)
  g_n_2 <- min(b_n_2 / d_n_2, 1)

  # compute the estimator
  estimate <- (1 - g_n_2) * cov_mat
  diag(estimate) <- diag(estimate) + g_n_2 * m_n
  estimate
}

linearShrinkLWEst(iris[1:4])
linearShrinkLWEst(iris[1:4] - 4)
#              Sepal.Length Sepal.Width Petal.Length Petal.Width
# Sepal.Length   0.68911527 -0.04211666     1.264785   0.5124098
# Sepal.Width   -0.04211666  0.19710838    -0.327191  -0.1207297
# Petal.Length   1.26478546 -0.32719104     3.101522   1.2859202
# Petal.Width    0.51240976 -0.12072969     1.285920   0.5852109

obj.bigsnp <- bigsnpr::snp_attach("tmp-data/GWAS_data_sorted_QC.rds")
G <- obj.bigsnp$genotypes
anyNA(dat <- G[, 1:500])

cov1 <- linearShrinkLWEst(dat)
cov2 <- linearShrinkLWEst(dat + 10)
all.equal(cov1, cov2)

cor1 <- cov2cor(cov1)
cor0 <- cor(dat)
plot(cor1, cor0)

microbenchmark::microbenchmark(
  PREV = matrixStats::sum2(cov_mat * diag(p)) / p,
  NEW = mean(diag(cov_mat)),
  check = "equal", times = 20
)
microbenchmark::microbenchmark(
  PREV = matrixStats::sum2((cov_mat - m_n * diag(p))^2) / p,
  NEW = norm(cov_mat, type = "F")^2 / p - m_n^2,
  check = "equal", times = 20
)

all.equal(cov(dat), crossprod(dat_centered) / (n - 1))

RSpectra::eigs_sym(cor0, k = 6)$values
# 16.275262 13.264047 12.589189 10.952114  9.733251  9.384796
RSpectra::eigs_sym(cor1, k = 6)$values
# 15.242677 12.292596 11.231468  9.801697  8.860538  8.703303
RSpectra::eigs_sym(cor0, k = 6, sigma = -0.5)$values
# 2.399382e-03  1.981330e-03  5.598204e-04  3.213755e-04  2.269980e-04 -1.443290e-15
RSpectra::eigs_sym(cor1, k = 6, sigma = -0.5)$values
# 0.05362454 0.05308275 0.05267710 0.05008856 0.04970662 0.04903219
