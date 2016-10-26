#' Date: 2016-10-24
#' Object: Test logistic regression
#' as iterated weighted linear regresssion
#' Results: OK, 17 min for all celiac

library(bigstatsr)
library(bigsnpr)
Rcpp::sourceCpp('tmp-tests/IRLS.cpp')

celiac2 <- AttachBigSNP(backingfile = "celiac_sub2_impute1",
                        backingpath = "../thesis-celiac/backingfiles")
bigX <- celiac2$genotypes
y <- (celiac2$fam$pheno + 1) / 2
X0 <- matrix(rnorm(2 * length(y)), ncol = 2)

X <- cbind(1, X0)
mod0 = glm(y ~ X - 1, family = binomial)
p0 = mod0$fitted
w0 = p0 * (1 - p0)
z0 = log(p0 / (1 - p0)) + (y - p0) / w0



N <- 20
X.sub <- sub.big.matrix(bigX, lastCol = N,
                        backingpath = celiac2$backingpath)

# 140 sec for step 0
print(system.time(
  test <- wcrossprod(X.sub@address, cbind(0, X), y, z0, w0,
                     tol = 1e-8, maxiter = 100)))

indNoConv <- which(!test$conv)
if ((l <- length(indNoConv)) > 0)
  printf(paste("For %d SNPs, IRLS has not converged",
               "using glm for those instead.\n", sep = "; "), l)
for (j in indNoConv) {
  mod <- glm(y ~ bigX[, j] + X - 1, family = binomial)
  coeffs <- summary(mod)$coefficients
  test$betas[j] <- coeffs[1]
  test$std[j] <- coeffs[2]
}

test2 <- matrix(0, 2, N)
for (i in 1:N) {
  mod <- glm(y ~ bigX[, i] + X - 1, family = binomial)
  test2[, i] = summary(mod)$coefficients[1, 1:2]
}
all.equal(test$betas, test2[1, ])
all.equal(test$std, test2[2, ])

# 17.5 min for all steps
print(system.time(
  test <- wcrossprod(bigX@address, cbind(0, X), y, z0, w0,
                     tol = 1e-8, maxiter = 100)))
