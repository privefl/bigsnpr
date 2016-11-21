#' Date: 2016-11-21
#' Object: Optimize pruning step
#' Results:

require(bigsnpr)

celiac3 <- AttachBigSNP("celiac-begin")
X <- celiac3$genotypes

require(bigstatsr)

size <- 50
thr <- 0.5

stats <- big_colstats(X)
keep <- rep(TRUE, ncol(X))
p <- stats$sum / (2 * nrow(X))
maf <- pmin(p, 1 - p)

test <- R_squared_chr2(X@address, keep,
                       maf, stats$sum, sqrt(stats$var),
                       size, thr)
ind <- which(test)

snps <- scan("../plink-1.07-x86_64/plink.prune.in", what = "character")
ind2 <- which(celiac3$map$marker.ID %in% snps)

all.equal(ind, ind2)
