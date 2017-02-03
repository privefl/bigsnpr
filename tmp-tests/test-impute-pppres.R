# path <- snp_readBed("../popres/POPRES_data/POPRES_allchr.bed",
#                     backingfile = "popres")
require(bigsnpr)
popresNA <- snp_attach("backingfiles/popresNA.bk", readonly = FALSE)

X <- popresNA$genotypes
X2 <- X[,]
indNA <- sort(sample(length(X), length(X) / 100))
X2[indNA] <- NA
X[] <- X2

print(system.time(
  test <- snp_impute(popresNA, ncores = 6, verbose = TRUE)
))
# 1h with ncores = 11
# 1h19 with ncores = 6

X3 <- test$genotypes

popres <- snp_attach("backingfiles/popres.bk")
mean(X3[indNA] != popres$genotypes[indNA])


# f <- MASS::fitdistr(nbNA, 'weibull')
# nbNA.simu <- round(rweibull(5e5, shape = 0.6, scale = 3))
