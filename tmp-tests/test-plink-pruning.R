#' Date: 2016-11-21
#' Object: Try to reproduce PLINK pruning procedure
#' Results: Seems OK, but will be very slooooooow..


require(bigsnpr)

# celiac <- AttachBigSNP("celiac_sub2_impute1", "../thesis-celiac/backingfiles/")
#
# celiac3 <- sub.bigSNP(celiac, ind.col = 1:5000)
#
# tmp <- BigToBed(celiac3, "tmp-files/celiac-begin.bed")
#
# celiac3 <- BedToBig("tmp-files/celiac-begin.bed",
#                     backingfile = "celiac-begin")

celiac3 <- AttachBigSNP("celiac-begin")
X <- celiac3$genotypes

# after plink --indep-pairwise 50 5 0.5 -> step 1 now
snps <- scan("../plink-1.07-x86_64/plink.prune.in", what = "character")
# deterministic -> OK

snps50 <- celiac3$map$marker.ID[1:50]
ind <- which(snps50 %in% snps)

X50 <- X[, 1:50]
corr <- cor(X50)
corr2 <- corr^2
ind2 <- which(corr2 > 0.5, arr.ind = TRUE)
ind3 <- ind2[ind2[, "row"] < ind2[, "col"], ]
ind.pruned <- unique(ind3[, "col"])

p <- colMeans(X50) / 2
maf <- pmin(p, 1 - p)
ind4 <- ifelse(maf[ind3[, "row"]] > maf[ind3[, "col"]], ind3[, "col"], ind3[, "row"])
ind5 <- sort(unique(ind4))
all.equal(ind5, setdiff(1:50, ind))


