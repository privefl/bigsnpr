Rcpp::sourceCpp('src/counts.cpp')
require(bigsnpr)

celiac <- AttachBigSNP("celiac_impute1_sub1")
#celiac <- Impute(celiac, ncores = 2, verbose = TRUE)
#celiac <- GetPhenos(celiac)
#test <- Counts(celiac)

X <- celiac$genotypes

indCase <- which(celiac$fam$pheno == 1)
indControl <- which(celiac$fam$pheno == -1)

#test <- mycount(X@address, indCase, indControl)

require(microbenchmark)

print(microbenchmark(
  test <- mycount(X@address, indCase, indControl),
  test2 <- mycount3(X@address, indCase, indControl),
  times = 20
))

# no significant difference
