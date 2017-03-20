require(bigsnpr)

test <- snp_attach("backingfiles/celiac300_sub1.rds")
G <- test$genotypes
G@code <- c(0, 1, rep(2, 254))

ind.train <- sort(sample(rows_along(G), size = 10e3))

ind.excl <- snp_indLRLDR(test$map$chromosome, test$map$physical.pos)
ind.keep <- snp_clumping(G, test$map$chromosome, thr.r2 = 0.2,
                         ind.row = ind.train, exclude = ind.excl)
svd <- big_randomSVD(G, snp_scaleBinom(), ind.col = ind.keep, ncores = 3)

plot(svd$u)
plot(svd$u[, 3:4])
plot(svd$d)

y01 <- bigstatsr:::transform_levels(test$fam$affection)
gwas <- big_univLogReg(G, y01.train = y01[ind.train], ind.train = ind.train,
                       covar.train = svd$u[ind.train, 1:2], ncores2 = 2)

manhattan2(gwas, test$map)
