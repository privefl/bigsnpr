snp <- snp_attachExtdata()
G <- snp$genotypes
obj.svd <- big_SVD(G, fun.scaling = snp_scaleBinom())
U <- obj.svd$u
mean(U[, 1, drop = FALSE])
all.equal(
  bigsnpr:::linRegPcadapt_cpp(G, U[, 1, drop = FALSE],
                              rows_along(G), cols_along(G))[, 1],
  big_univLinReg(G, U[, 1])$score
)
all.equal(
  bigsnpr:::linRegPcadapt_cpp(G, U[, 1, drop = FALSE],
                              rows_along(G), cols_along(G))[, 1],
  obj.svd$v[, 1] * obj.svd$d[1]
)


coef <- sapply(summary(lm(G[] ~ U[, 1])), function(mod) mod$coefficients[2, 3])
all.equal(coef, big_univLinReg(G, U[, 1])$score, check.attributes = FALSE,
          tolerance = 1e-18)

coef2 <- sapply(summary(lm(scale(G[]) ~ U[, 1])), function(mod) mod$coefficients[2, 3])
all.equal(coef2, big_univLinReg(G, U[, 1])$score, check.attributes = FALSE,
          tolerance = 1e-18)

coef3 <- sapply(summary(lm(scale(G[]) ~ U[, 1] + 0)), function(mod) mod$coefficients[3])
all.equal(coef3, big_univLinReg(G, U[, 1])$score, check.attributes = FALSE)

bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
obj.bed <- bed(bedfile)
test <- bigsnpr:::multLinReg(obj.bed, rows_along(G), cols_along(G),
                             U[, 1, drop = FALSE])
all.equal(test[, 1], big_univLinReg(G, U[, 1])$score, check.attributes = FALSE)

bed4pcadapt <- pcadapt::read.pcadapt(bedfile, type = "bed")
obj.pcadapt <- pcadapt::pcadapt(bed4pcadapt, K = 10, min.maf = 0)
test2 <- obj.pcadapt$zscores[, 1]
all.equal(test2, -big_univLinReg(G, U[, 1])$score)
plot(test2, -big_univLinReg(G, U[, 1])$score)
all.equal(obj.svd$v[, 1] * obj.svd$d[1], -test2)

ms <- big_scale()(G); ms$scale <- ms$scale * sqrt(nrow(G) - 1)
all.equal(ms$center, obj.svd$center)
test3 <- big_cprodMat(G, obj.svd$u, ncores = nb_cores(),
                      center = ms$center, scale = ms$scale)
all.equal(test3[, 1] * sqrt(nrow(G) - 2), big_univLinReg(G, U[, 1])$score)
test4 <- test3 / sqrt(1 - test3^2)
all.equal(test4[, 1] * sqrt(nrow(G) - 2), big_univLinReg(G, U[, 1])$score)


all.equal(obj.pcadapt$singular.values, obj.svd$d / sqrt((nrow(G) - 1) * ncol(G)))
