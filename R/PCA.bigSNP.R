
################################################################################



################################################################################

#'@name BigXXt
#'@title \eqn{X X^t} for genotype matrix of a "bigSNP".
#'@export
BigXXt <- function(x, block.size, ind.train = seq(nrow(x$genotypes)),
                   to.save = F) {
  n <- length(ind.train)
  if (to.save) {
    newfile <- checkFile(x, "PCA")
    bigK = bigmemory::big.matrix(n, n, type = "double", init = 0,
                                 backingfile = paste0(newfile, ".bk"),
                                 backingpath = x$backingpath,
                                 descriptorfile = paste0(newfile, ".desc"))
  } else {
    bigK = bigmemory::big.matrix(n, n, type = "double",
                                 init = 0, shared = F)
  }


  # compute p
  p.all <- bigcolsums((x$genotypes)@address, ind.train) / (2*n)

  # function to compute X*X^T
  BigXXt2 <- function(X) {
    printf("Computation of X * t(X)\n")
    intervals <- CutBySize(ncol(X), block.size)
    nb.block <- nrow(intervals)

    if (intr <- interactive()) {
      pb <- txtProgressBar(min = 0, max = nb.block, style = 3)
    }

    for (j in 1:nb.block) {
      if (intr) setTxtProgressBar(pb, j-1)
      ind <- seq2(intervals[j, ])
      incrSup(bigK@address, tcrossprod(center_p(X[ind.train, ind], p.all[ind])))
    }

    complete(bigK@address)

    if (intr) {
      setTxtProgressBar(pb, nb.block)
      close(pb)
    }

    return()
  }

  # compute K
  BigXXt2(x$genotypes)

  return(bigK)
}

################################################################################
