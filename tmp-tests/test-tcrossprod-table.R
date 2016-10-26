#' Date: 2016-10-26
#' Object: Test whether using a table for computations of tcrossprod
#' is faster?
#' Results: Not fast enough.


library(bigsnpr)
library(bigstatsr)

test <- AttachBigSNP("celiac_sub2_impute1", "../thesis-celiac/backingfiles")

X <- test$genotypes
means <- colmeans(X)
p <- means / 2
sds <- sqrt(2 * p * (1 - p))

tab2D <- sapply(0:2, function(x) (x - means) / sds)
tab3D <- apply(tab2D, 1, tcrossprod)
attributes(tab3D) <- list(dim = c(3, 3, ncol(tab3D)))



Rcpp::sourceCpp('tmp-tests/test-tcrossprod-table.cpp')
tabcrossprod(X@address, tab3D)

mat <- X[, 1:100]
print(typeof(mat)) # integer

n <- nrow(X)
res <- matrix(0, n, n) # 1 Gb

res <- tabcrossprod(mat, n, 100, res, tab3D)
