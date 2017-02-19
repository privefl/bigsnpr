require(bigsnpr)

popres <- snp_attach("backingfiles/popres.bk")
X <- popres$genotypes

prune.in <- scan("../../Téléchargements/plink_linux_x86_64/plink.prune.in",
                 what = "character")
test2 <- match(prune.in, popres$map$marker.ID)

print(system.time(
  test <- snp_pruning(popres, size = 100, is.size.in.kb = FALSE)
))

mean(test2 %in% test)
mean(test %in% test2)

ind <- which(!(test %in% test2))
ind2 <- test[ind]

ind3 <- which(!(test2 %in% test))
ind4 <- test2[ind3]

X.cor <- cor(X[, 1:100])
ms <- colMeans(popres$genotypes[, ind2])
sds <- ms <- apply(popres$genotypes[, ind2], 2, sd)

dist.chr1 <- popres$map$physical.pos[popres$map$chromosome == 1]
dist.size50 <- diff(dist.chr1, lag = 50) / 1000
dist.size100 <- diff(dist.chr1, lag = 100) / 1000

head(test)
test2 <- match(prune.in, popres$map$marker.ID)
head(test2)
mean(test2 %in% test)

plot(dist.chr1)
size <- 250e3
first <- min(which(dist.chr1 - dist.chr1[1] > size))
dist.chr1 <-
