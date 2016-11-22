#' Date: 2016-11-21
#' Object: Optimize pruning step
#' Results:

require(bigsnpr)

celiac2 <- AttachBigSNP("celiac_sub2_impute1",
                        "../thesis-celiac/backingfiles/")
print(dim(celiac2$genotypes))

# print(system.time(BigToBed(celiac2, "tmp-files/celiac-imputed.bed")))
# 50 sec

print(system.time(
  test <- PrunePlink(celiac2, size = 50, thr.corr = 0.5, ncores = 6)
)) # 38 sec versus 2h30

# after plink --indep-pairwise 50 5 0.5
snps <- scan("../plink-1.07-x86_64/test-celiac.prune.in", what = "character")
ind2 <- which(celiac2$map$marker.ID %in% snps)
print(mean(ind2 %in% test)) # 99.6%
