#' Date: 2016-11-21
#' Object: Optimize pruning step
#' Results:

require(bigsnpr)

celiac3 <- AttachBigSNP("celiac-begin")
X <- celiac3$genotypes
print(dim(X))

require(bigstatsr)
Rcpp::sourceCpp('tmp-tests/test-plink-pruning.cpp')

prune2 <- function(X, size = 50, thr = 0.5) {
  stats <- big_colstats(X)
  keep <- rep(TRUE, ncol(X))
  n <- nrow(X)
  p <- stats$sum / (2 * n)
  maf <- pmin(p, 1 - p)
  denoX <- (n - 1) * stats$var

  keep <- R_squared_chr2(X@address,
                         keep = keep,
                         mafX = maf,
                         sumX = stats$sum,
                         denoX = denoX,
                         size = size,
                         thr = thr)
}

print(system.time(test <- prune2(X))) # 2 sec versus 8 min
ind <- which(test)

# after plink --indep-pairwise 50 1 0.5
snps <- scan("../plink-1.07-x86_64/plink.prune.in", what = "character")
ind2 <- which(celiac3$map$marker.ID %in% snps)
print(all.equal(ind, ind2))
