#' Date: 2016-11-21
#' Object: Try to reproduce PLINK pruning procedure
#' Results:


require(bigsnpr)

# celiac <- AttachBigSNP("celiac_sub2_impute1", "../thesis-celiac/backingfiles/")
#
# celiac3 <- sub.bigSNP(celiac, ind.col = 1:5000)
#
# tmp <- BigToBed(celiac3, "celiac-begin.bed")
#
# write(celiac3$map$marker.ID, "snps.txt")

celiac3 <- BedToBig("tmp-files/celiac-plink.bed",
                    backingfile = "celiac-begin",
                    block.size = 1000)
X <- celiac3$genotypes

# after plink --indep-pairwise 50 5 0.5
snps <- scan("../plink-1.07-x86_64/plink.prune.in", what = "character")

snps50 <- celiac3$map$marker.ID[1:50]
ind <- which(snps50 %in% snps)

X50 <- X[, 1:50]
corr <- cor(X50, use = "pairwise.complete.obs")
corr2 <- corr^2
ind2 <- which(corr2 > 0.5, arr.ind = TRUE)
ind3 <- ind2[ind2[, "row"] < ind2[, "col"], ]
ind.pruned <- unique(ind3[, "col"]) # not that..
