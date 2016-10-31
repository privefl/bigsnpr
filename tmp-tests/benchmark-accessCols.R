Rcpp::sourceCpp('src/counts.cpp')
require(bigsnpr)

celiac <- AttachBigSNP("test")
celiac <- GetPhenos(celiac)
test <- Counts(celiac)

X <- celiac$genotypes

indCase <- sort(sample(nrow(X), 5000))
indControl <- setdiff(seq(nrow(X)), indCase)

require(microbenchmark)

print(microbenchmark(
  test <- mycount(X@address, indCase, indControl),
  test2 <- mycount3(X@address, indCase, indControl),
  times = 10
))

