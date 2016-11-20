#' Date: 2016-10-25
#' Object: Testing GRM two-step for relationship
#' Results:

library(bigsnpr)

celiac2 <- AttachBigSNP(backingfile = "celiac_impute1_sub1")
X <- celiac2$genotypes

nb <- 5000
ind <- round(1:nb / (nb + 1) * ncol(X))
X2 <- deepcopy(X, rows = 1:2000, cols = ind, shared = FALSE)

require(bigstatsr)

f1 <- function(X, ind.train) {
  means <- big_colstats(X, ind.train)$sum / length(ind.train)
  p <- means / 2
  sds <- sqrt(2 * p * (1 - p))
  list(mean = means, sd = sds)
}
print(system.time(test <- big_tcrossprodSelf(X2, f1)))

ind <- which(test[,] / ncol(X2) > 0.05, arr.ind = TRUE)
ind2 <- ind[ind[, "row"] > ind[, "col"], ]
ind3 <- sort(unique(as.vector(ind2)))

print(system.time(test2 <- big_tcrossprodSelf(X, f1, ind.train = ind3)))

ind <- which(test2[,] / ncol(X) > 0.05, arr.ind = TRUE)
ind2 <- ind[ind[, "row"] > ind[, "col"], ]
print(test2[ind2] / ncol(X))
