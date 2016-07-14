#'@title Quality control, subset and imputation for a "bigSNP".
#'@name qcimpute
#'@param x A \code{bigSNP}.
#'@return A new \code{bigSNP}.
#'@examples #TODO
NULL

################################################################################

#'@name QC
#'@description \code{QC}: Quality control (filters)
#'for a \code{bigSNP} resulting
#'in a \code{bigSNP} of lower dimension.
#'@rdname qcimpute
#'@seealso \code{\link{bigSNP-class}}
#'@export
QC <- function(x) {
  if (class(x) != "bigSNP") stop("x must be a bigSNP")

  printf("Not yet implemented\n")
}

#'@description \code{sub.bigSNP}: a function
#'to get a subset of an object of class \code{bigSNP}.
#'@param ind.row Indices of the rows (individuals) to keep.
#'@param ind.col Indices of the columns (SNPs) to keep.
#'@export
#'@name sub.bigSNP
#'@rdname qcimpute
sub.bigSNP <- function(x, ind.row, ind.col) {
  if (class(x) != "bigSNP") stop("x must be a bigSNP")

  number <- 1
  while (file.exists(
    file.path(x$backingpath,
              paste0(newfile <- paste0(x$backingfile, "_sub", number),
                     ".desc")))) {
    number <- number + 1
  }
  X2 <- bigmemory::big.matrix(length(ind.row), length(ind.col), type = "char",
                              backingfile = paste0(newfile, ".bk"),
                              backingpath = x$backingpath,
                              descriptorfile = paste0(newfile, ".desc"))

  deepcopyPart((x$genotypes)@address, X2@address, ind.row, ind.col)

  snp_list <- list(genotypes = X2,
                   fam = x$fam[ind.row, ],
                   map = x$map[ind.col, ],
                   backingfile = newfile,
                   backingpath = x$backingpath)
  class(snp_list) <- "bigSNP"

  saveRDS(snp_list, file.path(x$backingpath, paste0(newfile, ".rds")))

  return(snp_list)
}

################################################################################

#'@description \code{Impute}: Imputation function
#'for a \code{bigSNP}.
#'@export
#'@name Impute
#'@param ncores Number or cores used.
#'Default doesn't use parallelism.
#'@rdname qcimpute
Impute <- function(x, ncores = 1) {
  if (class(x) != "bigSNP") stop("x must be a bigSNP")

  # get descriptors
  X.desc <- describe(x$genotypes)

  number <- 1
  while (file.exists(
    file.path(x$backingpath,
              paste0(newfile <- paste0(x$backingfile, "_impute", number),
                     ".desc")))) {
    number <- number + 1
  }
  X2 <- deepcopy(x$genotypes, type = "char",
                 backingfile = paste0(newfile, ".bk"),
                 backingpath = x$backingpath,
                 descriptorfile = paste0(newfile, ".desc"))
  X2.desc <- describe(X2)

  # function that imputes one chromosome
  ImputeChr <- function(lims) {
    #printf("Imputing chromosome %d with \"nearest neighbors\"...\n", lims[3])

    X <- sub.big.matrix(X.desc, firstCol = lims[1], lastCol = lims[2],
                             backingpath = x$backingpath)
    X2 <- sub.big.matrix(X2.desc, firstCol = lims[1], lastCol = lims[2],
                             backingpath = x$backingpath)

    predictNA <- function(ind, ind2) {
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

          X2[i, ind[ind2]] <- round(pred)
        }
      }

      return(0)
    }

    m <- ncol(X)

    opt.save <- options(bigmemory.typecast.warning = FALSE)

    # first three columns
    for (j in 1:3) {
      predictNA(1:7, j)
    }

    # middle
    for (j in 4:(m-3)) {
      predictNA(j + -3:3, 4)
    }

    # last three columns
    for (j in 5:7) {
      predictNA(m + -6:0, j)
    }

    options(opt.save)

    return(0)
  }

  range.chr <- LimsChr(x)

  obj <- foreach::foreach(i = 1:nrow(range.chr),
                          .noexport = c("x", "X2"))
  expr_fun <- function(i) ImputeChr(range.chr[i, ])
  foreach2(obj, expr_fun, ncores)

  snp_list <- list(genotypes = X2,
                   fam = x$fam,
                   map = x$map,
                   backingfile = newfile,
                   backingpath = x$backingpath)
  class(snp_list) <- "bigSNP"

  saveRDS(snp_list, file.path(x$backingpath, paste0(newfile, ".rds")))

  return(snp_list)
}

################################################################################
