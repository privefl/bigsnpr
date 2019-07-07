snp <- snp_attachExtdata()
G <- snp$genotypes
obj.svd <- big_SVD(G, fun.scaling = snp_scaleBinom())
U <- obj.svd$u
mean(U[, 1, drop = FALSE])
plot(
  bigsnpr:::multLinReg(G, rows_along(G), cols_along(G), U[, 1, drop = FALSE])[, 1],
  big_univLinReg(G, U[, 1])$score
)
all.equal(
  bigsnpr:::multLinReg(G, rows_along(G), cols_along(G), U),
  sapply(1:10, function(k) big_univLinReg(G, U[, k])$score)
)
system.time(sapply(1:10, function(k) big_univLinReg(G, U[, k])$score)) /
  system.time(bigsnpr:::multLinReg(G, rows_along(G), cols_along(G), U)) # 10-12

# coef <- sapply(summary(lm(G[] ~ U[, 1])), function(mod) mod$coefficients[2, 3])
# all.equal(coef, big_univLinReg(G, U[, 1])$score, check.attributes = FALSE,
#           tolerance = 1e-18)
#
# coef2 <- sapply(summary(lm(scale(G[]) ~ U[, 1])), function(mod) mod$coefficients[2, 3])
# all.equal(coef2, big_univLinReg(G, U[, 1])$score, check.attributes = FALSE,
#           tolerance = 1e-18)
#
# coef3 <- sapply(summary(lm(scale(G[]) ~ U[, 1] + 0)), function(mod) mod$coefficients[3])
# all.equal(coef3, big_univLinReg(G, U[, 1])$score, check.attributes = FALSE)

bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
bed4pcadapt <- pcadapt::read.pcadapt(bedfile, type = "bed")
obj.pcadapt <- pcadapt::pcadapt(bed4pcadapt, K = 10, min.maf = 0)
# all.equal(obj.pcadapt$singular.values, obj.svd$d / sqrt((nrow(G) - 1) * ncol(G)))
test2 <- obj.pcadapt$zscores[, 1]
all.equal(test2, -big_univLinReg(G, U[, 1])$score)
plot(test2, -big_univLinReg(G, U[, 1])$score)
all.equal(obj.svd$v[, 1] * obj.svd$d[1], -test2)
all.equal(obj.svd$v[, 1] * obj.svd$d[1], big_univLinReg(G, U[, 1])$score)

ms <- big_scale()(G); ms$scale <- ms$scale * sqrt(nrow(G) - 1)
all.equal(ms$center, obj.svd$center)
test3 <- big_cprodMat(G, obj.svd$u, ncores = nb_cores(),
                      center = ms$center, scale = ms$scale)
all.equal(test3[, 1] * sqrt(nrow(G) - 2), -test2)
all.equal(test3[, 1] * sqrt(nrow(G) - 2), big_univLinReg(G, U[, 1])$score)
test4 <- test3 / sqrt(1 - test3^2)
all.equal(test4[, 1] * sqrt(nrow(G) - 2), big_univLinReg(G, U[, 1])$score)
