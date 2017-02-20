# Convert bed/bim/fam to vcf
# ./plink -recode vcf -bfile ../../Bureau/POPRES_data/popresSub

# Impute with BEAGLE 4.1
# java -Xmx6g -jar ../beagle.21Jan17.6cc.jar gt=plink.vcf out=test


require(bigsnpr)

popres <- snp_attach("backingfiles/popres.bk")
X <- popres$genotypes

maf <- snp_MAF(X)
ind <- which(popres$map$chromosome == 1 & maf > 0.05)

popres.chr1 <- sub.bigSNP(popres, ind.col = ind)
X2 <- popres.chr1$genotypes
n <- nrow(X2)
m <- ncol(X2)

print(table(
  nbNA <- VGAM::rbetabinom.ab(m, size = n, shape1 = 0.6, shape2 = 50)
))

for (j in 1:m) {
  indNA <- sample(n, size = nbNA[j])
  X2[indNA, j] <- NA
}
stopifnot(all.equal(colSums(is.na(X2[,])), nbNA))

print(system.time(
  test <- snp_impute(popres.chr1)
))
# 15 min

popres.chr1.noNA <- snp_attach("backingfiles/popres_sub1.bk")
indNA <- which(is.na(X2[,]))
mean(test$genotypes[indNA] != popres.chr1.noNA$genotypes[indNA]) # error: 4.26%
infos <- test$imputation
print(sum(infos$nbNA * infos$error, na.rm = TRUE) / length(indNA)) # estimated: 4.23%
plot(infos, pch = 19, cex = 0.5)
curve(5/x, col = "blue", add = TRUE)
curve(10/x, col = "red", add = TRUE)


# beagle time
snp_writeBed(popres.chr1, "../../plink_linux_x86_64/popres_chr1_NA.bed")
