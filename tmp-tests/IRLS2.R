#' Date: 2016-10-25
#' Object: Testing GWAS with ncores = 10
#' Results: 3.5 min :D

library(bigsnpr)

celiac2 <- AttachBigSNP(backingfile = "celiac_sub2_impute1",
                        backingpath = "../thesis-celiac/backingfiles")

n <- nrow(celiac2$genotypes)
X0 <- matrix(rnorm(2 * n), ncol = 2)

print(system.time(
  test <- GWAS(celiac2, covar = X0, ncores = 10)))

pval <- 2 * pnorm(-abs(test[1, ] / test[2, ]))

lpval <- -log10(pval)
plot(replace(lpval, lpval < 1, NA), type = "h", col = celiac2$map$chromosome)
