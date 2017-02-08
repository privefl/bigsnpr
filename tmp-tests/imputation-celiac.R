require(bigsnpr)

celiac <- snp_attach("../thesis-celiac/backingfiles/celiac.bk")

print(system.time(
  test <- snp_impute(celiac, ncores = 11, verbose = TRUE)
)) # 112 min

info <- test$imputation
plot(info, cex = 0.5, pch = 19, xlim = c(0, 600))
curve(10/x, col = "red", add = TRUE, from = 1e-6, to = max(info$nbNA), n = 1e4)
curve(20/x, col = "blue", add = TRUE, from = 1e-6, to = max(info$nbNA), n = 1e4)

ind <- which(info$nbNA > 500 & info$error> 0.4)
ind2 <- which(info$nbNA * info$error > 20 | info$nbNA > 600)
ind3 <- which(info$nbNA * info$error > 12 | info$nbNA > 600)

X2 <- test$genotypes[, info$error > 0.6]
colMeans(X2)

sum(is.na(test$genotypes[,]))

ind4 <- which(info$nbNA == 0)
ind5 <- union(ind2, ind4)
sum(info[-ind5, "nbNA"] * info[-ind5, "error"]) / sum(info[-ind2, "nbNA"])
