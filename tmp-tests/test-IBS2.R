require(bigsnpr)
Rcpp::sourceCpp('src/IBS.cpp')

test <- snp_attach("backingfiles/popresNA_sub1.rds")
G <- test$genotypes
n <- nrow(G)

resIBS <- tmpFBM(4 * n, n, type = "integer", init = 0)

table <- cbind(rbind(matrix(c(2, 1, 0, 1, 2, 1, 0, 1, 2), 3), 3), 3)

print(system.time(
  IBS(attach.big.matrix(G)@address,
      code = c(0, 1, 2, rep(3, 253)),
      table = table,
      attach.big.matrix(resIBS)@address)
))

attach.BM(resIBS)[1:5, 1:5]

resIBS2 <- array(data = 0L, dim = c(4, n, n))

print(system.time(
  resIBS2 <- IBS2(attach.big.matrix(G)@address,
                  code = c(0, 1, 2, rep(3, 253)),
                  table = table,
                  resIBS2)
))

G@code <- bigsnpr:::CODE_DOSAGE
print(system.time(
  xxt <- big_tcrossprodSelf(G, snp_scaleBinom())
))
