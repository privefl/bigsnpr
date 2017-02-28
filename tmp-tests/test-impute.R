require(bigsnpr)

popres.chr1 <- snp_attach("backingfiles/popres_sub2.bk")
X2 <- popres.chr1$genotypes
n <- nrow(X2)
m <- ncol(X2)


print(system.time(
  test <- snp_imputeCV2(popres.chr1)
))
# 2h12 with default parameters -> 4.03%

# 18h with cross-validation with high-parameters -> 4.21%
# snp_imputeCV(popres.chr1, nrounds = 50, max_depth = 6,
#              sizes = c(20, 30, 50, 100), Kfolds = c(5, 6, 8, 10))

# 2h40 -> 4.3%
# snp_imputeCV(popres.chr1, nrounds = 10, max_depth = 3,
#              sizes = c(20, 30, 50, 100), Kfolds = c(5, 6, 8, 10)

# 1h20 -> 4.4%
# test <- snp_imputeCV(popres.chr1, nrounds = 10, max_depth = 3,
#                      sizes = c(20, 30, 50, 100), Kfolds = rep(5, 4))

# 46 min -> 4.3%
# snp_imputeCV(popres.chr1, nrounds = 10, max_depth = 3, Kfolds = rep(5, 4))

# 1h40 -> 4.03%
# snp_imputeCV(popres.chr1)

popres.chr1.noNA <- snp_attach("backingfiles/popres_sub1.bk")
indNA <- which(is.na(X2[,]))
mean(test$genotypes[indNA] != popres.chr1.noNA$genotypes[indNA]) # error: 4.03%
infos <- test$imputation
print(sum(infos$nbNA * infos$error, na.rm = TRUE) / length(indNA)) # estimated: 4.14%
plot(infos, pch = 19, cex = 0.5)
curve(5/x, col = "blue", add = TRUE)
curve(10/x, col = "red", add = TRUE)

require(RcppRoll)
nbNA_roll <- roll_mean(infos$nbNA, n = 20, fill = NA)
plot(infos$nbNA, nbNA_roll)
plot(nbNA_roll, infos$error)
abline(lm(infos$error ~ nbNA_roll), col = "red")
ind <- which(infos$error > 0.3)

