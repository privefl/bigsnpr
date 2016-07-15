
################################################################################



################################################################################

#'@name BigXXt
#'@description Compute \eqn{X X^T} for the genotype matrix of a "bigSNP"
#'with a particular scaling.
#'@title \code{Tcrossprod} for the genotype matrix of a "bigSNP".
#'@param x A \code{bigSNP}.
#'@param block.size Maximum number of loci read at once (for all individuals).
#'@param ind.train An optional vector of the row indices that are used.
#' If not specified, all data are used.
#'@param to.save Is the result filebacked ? Default is \code{FALSE}.
#'@seealso \code{\link{bigSNP}}
#'@return A \code{big.matrix} of type \code{double}.
#'@export
BigXXt <- function(x, block.size, ind.train = seq(nrow(x$genotypes)),
                   to.save = FALSE) {
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
  p.all <- bigcolsumsChar((x$genotypes)@address, ind.train) / (2*n)

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
      p.ind <- p.all[ind]
      mean <- 2*p.ind
      sd <- sqrt(2*p.ind*(1-p.ind))
      incrSup(bigK@address, tcrossprod(
        scaling(X[ind.train, ind], mean, sd)))
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

PCA.bigSNP <- function(x, block.size, k = NULL,
                       ind.train = seq(nrow(x$genotypes)),
                       to.save = FALSE) {
  bigK <- BigXXt(x, block.size, ind.train, to.save)
  n <- nrow(bigK)

  means <- bigcolsumsDouble(bigK, ind.train) / n
  symCenter(K2, means, mean(means))

  if (is.null(k)) {
    eig <- eigen(bigK[,], symmetric = TRUE)
  } else {
    eig <- RSpectra::eigs_sym(bigK[,], k)
  }

  alphas <- scaling(eig$vectors,
                    rep(0, length(eig$values)),
                    sqrt(eig$values))
  rm(eig)
}

################################################################################
