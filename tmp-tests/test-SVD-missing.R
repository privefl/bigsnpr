library(bigsnpr)
obj.bed <- bed("tmp-data/1000G_phase3_common_hapmap_norel.bed")
set.seed(1); ind.col <- seq(1, ncol(obj.bed), by = 20)
obj.svd <- bed_randomSVD(obj.bed, ind.col = ind.col)
U <- obj.svd$u

# snp_readBed2("tmp-data/1000G_phase3_common_hapmap_norel.bed", "tmp-data/1000G",
#              ind.col = ind.col)
snp <- snp_attach("tmp-data/1000G.rds")
G <- snp$genotypes

set.seed(NULL)
G2 <- big_copy(G); M <- ncol(G); N <- nrow(G)
hist(nbNA <- VGAM::rbetabinom.ab(M, size = N, shape1 = 0.7, shape2 = 10))
sum(nbNA) / length(G)
indNA <- cbind(
  unlist(
    lapply(nbNA, function(nb) {
      `if`(nb > 0, sample(N, size = nb), NULL)
    })
  ),
  rep(cols_along(G), nbNA)
)
G2[indNA] <- 3
round(100 * big_counts(G2)[4, ] / nrow(G2), 1)
G2[, 1]
snp2 <- snp
snp2$genotypes <- G2
bedfile2 <- snp_writeBed(snp2, tempfile(fileext = ".bed"))

obj.bed2 <- bed(bedfile2)
svd2 <- bed_randomSVD(obj.bed2)
U2 <- svd2$u
plot(U, U2)
plot(U[, 1:2], U2[, 1:2], pch = 20); abline(0, 1, col = "red")

svd3 <- bed_randomSVD(obj.bed2, ncores = 2)
U3 <- svd3$u
plot(U, U3)
plot(U[, 1:2], U3[, 1:2], pch = 20); abline(0, 1, col = "red")

all.equal(U2, U3) # basically the same
round(100 * cor(U2, U), 1)
round(100 * cor(U3, U), 1)

all.equal(obj.svd$d, svd2$d)
all.equal(obj.svd$d, svd3$d)
obj.svd$d / svd2$d

K <- bed_tcrossprodSelf(obj.bed, ind.col = ind.col)
val <- eigen(K[], symmetric = TRUE, only.values = TRUE)$values
rbind(obj.svd$d, sqrt(head(val, length(obj.svd$d))))

K2 <- bed_tcrossprodSelf(obj.bed2)
val2 <- eigen(K2[], symmetric = TRUE, only.values = TRUE)$values
rbind(svd2$d, sqrt(head(val2, length(svd2$d))))
