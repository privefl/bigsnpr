snp_tcrossprodSelf <- function(X.,
                               ind.row = rows_along(X.),
                               ind.col = cols_along(X.),
                               block.size = 1000) {
  X <- attach.big.matrix(X.)
  n <- length(ind.row)
  K12 <- matrix(0, n, n)
  K34 <- matrix(0, n, n)
  K5 <- matrix(0, n, n)

  intervals <- bigstatsr:::CutBySize(length(ind.col), block.size)
  nb.block <- nrow(intervals)

  code1 <- c(1, 0, 0, rep(0, 253))
  code2 <- c(0, 0, 1, rep(0, 253))
  code3 <- c(1, 0, 1, rep(0, 253))
  code4 <- c(0, 1, 0, rep(0, 253))
  code5 <- c(1, -1, 1, rep(0, 253))

  for (j in 1:nb.block) {
    X.raw <- X[ind.row, ind.col[bigstatsr:::seq2(intervals[j, ])]]

    X1 <- decode(X.raw, code1)
    X2 <- decode(X.raw, code2)
    K12 <- bigstatsr:::incrMat(K12, tcrossprod(X1, X2))

    X3 <- decode(X.raw, code3)
    X4 <- decode(X.raw, code4)
    K34 <- bigstatsr:::incrMat(K34, tcrossprod(X3, X4))

    X5 <- decode(X.raw, code5)
    K5 <- bigstatsr:::incrSup2(K5, tcrossprod(X5))
  }

  K
}

decode <- bigstatsr:::decodeMat
X. <- G
test <- snp_tcrossprodSelf(G, snp_scaleBinom())



require(SNPRelate)

genofile <- snpgdsOpen("../../Bureau/POPRES_data/popresSub.gds")

print(system.time(
  relate <- snpgdsIBSNum(genofile, num.thread = 1,
                         snp.id = test$map$marker.ID)
)) # 7.5 min (6 cores) -> 12 min with PLINK (11 cores)
K0 <- K12 + t(K12)
K1 <- K34 + t(K34)
K5 <- bigstatsr:::complete2(K5)
K2 <- K5 + K1 - K0
all.equal(K0, relate$ibs0)
all.equal(K1, relate$ibs1)
all.equal(K2, relate$ibs2)

