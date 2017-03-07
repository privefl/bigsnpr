require(bigsnpr)

celiac <- snp_attach("backingfiles/testRead.rds")

print(system.time(
  test <- big_counts(celiac$genotypes)
)) # 38 sec -> 9 sec the second time

dim(test)

system.time(colSums(test))

print(system.time(
  test2 <- big_counts(celiac$genotypes, byrow = TRUE)
)) # 9 sec

dim(test2)

tmp <- as.numeric(rownames(test2))
iNA <- which(is.na(tmp))
test2[iNA, ]
rowSums(test2)
rowSums(test)

print(system.time(
  test4 <- snp_MAX3(celiac)
))

m <- length(test4$pS)
plot(-log10((m:1)/(m+1)), sort(-log10(test4$pS)), pch = 19, cex = 0.5,
     ylim = c(0, 20))
abline(0, 1, col = "red")

print(system.time(
  test5 <- snp_MAX3(celiac, val = 0.5)
))

m <- length(test5$pS)
plot(-log10((m:1)/(m+1)), sort(-log10(test5$pS)), pch = 19, cex = 0.5,
     ylim = c(0, 20))
abline(0, 1, col = "red")
