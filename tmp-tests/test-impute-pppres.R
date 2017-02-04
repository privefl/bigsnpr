pacman::p_load_current_gh("privefl/bigsnpr")

path2popres <- "../POPRES_data/POPRES_allchr.bed"

pathNA <- snp_readBed(path2popres, backingfile = "popresNA")
popresNA <- snp_attach(pathNA, readonly = FALSE)

X <- popresNA$genotypes
n <- nrow(X)
m <- ncol(X)

print(table(
  nbNA <- rbetabinom.ab(m, size = n, shape1 = 0.6, shape2 = 100)
))

for (j in 1:m) {
  indNA <- sample(n, size = nbNA[j])
  X[indNA, j] <- NA
}

print(system.time(
  test <- snp_impute(popresNA, ncores = 3, verbose = TRUE)
))
# 2h20 with 3 cores on laptop


pathNoNA <- snp_readBed(path2popres, backingfile = "popres")
popresNoNA <- snp_attach(pathNoNA)
X2 <- popresNoNA$genotypes

popresImpute <- snp_attach("backingfiles/popresNA_impute1.bk")
X3 <- popresImpute$genotypes

indNA <- which(is.na(X[,]), arr.ind = TRUE)
stopifnot(nrow(indNA) == sum(nbNA))
print(mean(X2[indNA] != X3[indNA])) # 3.4 %
info <- popresImpute$imputation
print(crossprod(info$nbNA, info$error) / sum(nbNA)) # 3.4 %
