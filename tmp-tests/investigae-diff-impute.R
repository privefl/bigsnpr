require(bigsnpr)

popres.XGB <- snp_attach("backingfiles/popres_sub2_impute1.bk")
popres.BEAGLE <- snp_attach("backingfiles/popres_impute_beagle.bk")
popres.NA <- snp_attach("backingfiles/popres_sub2.bk")
popres.NoNA <- snp_attach("backingfiles/popres_sub1.bk")

X.NA <- popres.NA$genotypes
m <- ncol(X.NA)
indNA <- which(is.na(X.NA[,]))
nbNA <- colSums(is.na(X.NA[,]))
indNA.arr <- which(is.na(X.NA[,]), arr.ind = TRUE)

XGB.true <- (popres.XGB$genotypes[indNA] == popres.NoNA$genotypes[indNA])
BEAGLE.true <- (popres.BEAGLE$genotypes[indNA] == popres.NoNA$genotypes[indNA])
sum(XGB.true)
sum(BEAGLE.true)
length(indNA)

ind.diff <- BEAGLE.true - XGB.true
cols.diff <- indNA.arr[ind.diff, "col"]
tab.diff <- table(cols.diff)

require(dplyr)
df_test <- data.frame(cbind(indNA.arr, ind.diff))
diff.col <- df_test %>% group_by(col) %>% summarise(Sum = sum(ind.diff))

nbNA.BeagleBest <- nbNA[diff.col$col[which(diff.col$Sum > 1)]]
nbNA.XGBBest <- nbNA[diff.col$col[which(diff.col$Sum < -1)]]

par.save <- par(mfrow = c(2, 1))
# hist(nbNA, freq = FALSE, xlim = c(0, 200), breaks = 100)
hist(nbNA.BeagleBest, freq = FALSE, xlim = c(0, 150), breaks = 100)
hist(nbNA.XGBBest, freq = FALSE, xlim = c(0, 150), breaks = 20)
par(par.save)


max.NA <- max(nbNA.BeagleBest, nbNA.XGBBest)
tab.nbNA.BeagleBest <- table(nbNA.BeagleBest)
tab.nbNA.XGBBest <- table(nbNA.XGBBest)
vec.nbNA.BeagleBest <- numeric(max.NA)
vec.nbNA.BeagleBest[as.numeric(names(tab.nbNA.BeagleBest))] <- as.vector(tab.nbNA.BeagleBest)
vec.nbNA.XGBBest <- numeric(max.NA)
vec.nbNA.XGBBest[as.numeric(names(tab.nbNA.XGBBest))] <- as.vector(tab.nbNA.XGBBest)

vecDiffratio <- log((vec.nbNA.BeagleBest + 0.5) / (vec.nbNA.XGBBest + 0.5))
plot(vecDiffratio)
