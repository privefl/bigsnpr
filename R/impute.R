#'@title Quality control, subset and imputation for a "bigSNP".
#'@name qcimpute
#'@param bigsnp A \code{bigSNP}.
#'@param backingfile The root name for the backing file(s) for the cache of
#'the resulting object.
#'@param backingpath The path to the directory containing the file backing cache.
#'Default is "backingfiles". It needs to exist.
#'@return A \code{bigSNP}.
#'@examples #TODO
NULL

################################################################################

#'@name QC
#'@description \code{QC}: Quality control (filters)
#'for a \code{\link{bigSNP}} resulting
#'in a \code{bigSNP} of lower dimension.
#'@rdname qcimpute
#'@export
QC <- function(bigsnp) {
  if (class(bigsnp) != "bigSNP") stop("bigsnp must be a bigSNP")
}

#'@description \code{subset.bigSNP}: an S3 generic function
#'for \code{subset} with the class \code{bigSNP}.
#'@param ind.row Indices of the rows (individuals) to keep.
#'@param ind.col Indices of the columns (SNPs) to keep.
#'@export
#'@method subset bigSNP
#'@name subset.bigSNP
#'@rdname qcimpute
subset.bigSNP <- function(bigsnp, ind.row, ind.col,
                          backingfile,
                          backingpath = "backingfiles") {

}

################################################################################

impute <- function(X.desc, lims) {
  bigMat <- sub.big.matrix(X.desc, firstCol = lims[1], lastCol = lims[2],
                           backingpath = "backingfiles")

  predictNA <- function(X, ind, ind2) {
    tmp <- X[, ind]
    indNA <- which(is.na(tmp[, ind2]))
    if (length(indNA) > 0) {
      for (i in indNA) {
        tmpSum <- rowSums(sweep(tmp[, -ind2], 2, tmp[i, -ind2], '=='), na.rm = T)

        k <- 6
        cond <- T
        while (cond) {
          indNN <- which(tmpSum == k)
          k <- k - 1
          pred <- mean(tmp[indNN, ind2], na.rm = T)
          cond <- is.na(pred)
        }

        X[i, ind[ind2]] <- round(pred)
      }
    }

    return(0)
  }

  m = ncol(bigMat)

  opt.save = options(bigmemory.typecast.warning = FALSE)

  # first three columns
  for (j in 1:3) {
    predictNA(bigMat, 1:7, j)
  }

  # middle
  for (j in 4:(m-3)) {
    predictNA(bigMat, j + -3:3, 4)
  }

  # last three columns
  near = bigMat[, m + -6:0]
  for (j in 5:7) {
    predictNA(bigMat, m + -6:0, j)
  }

  options(opt.save)

  return()
}

#'@description \code{impute}: Imputation function
#'for a \code{bigSNP}.
#'@export
#'@name impute
#'@rdname qcimpute
impute_all <- function(X, celiac, ncores) {
  range.chr <- foreach(chr = 1:22, .combine = 'rbind') %do% {
    range(which(celiac$map$chromosome[-celiac$indQC] == chr))
  }
  X.desc = describe(X)

  is.seq <- (ncores == 1)
  if (is.seq) {
    registerDoSEQ()
  } else {
    cl <- makeCluster(ncores, outfile = "")
    registerDoParallel(cl)
  }
  foreach(chr = 1:nrow(range.chr),
          .combine = 'c',
          .packages = c("bigmemory"),
          .export = c("impute", "printf")) %dopar% {
            printf("Imputing chromosome %d with nearest neighbors...\n", chr)
            impute(X.desc, range.chr[chr, ])
          }
  if (!is.seq) stopCluster(cl)

  return()
}
