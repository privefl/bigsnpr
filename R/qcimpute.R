#' @title Imputation, quality control and subset for a "bigSNP".
#' @name impute-qc-sub
#' @inheritParams bigsnpr-package
#' @return A new \code{bigSNP}.
#' @seealso \code{\link{bigSNP}}
#' @example examples/example.sub.bigSNP.R
NULL


################################################################################

#' @name QC
#' @description \code{QC}: Quality control (filters)
#' for a \code{bigSNP} resulting
#' in a \code{bigSNP} of lower dimension.
#' @param hwe.pval Level threshold (allowed type-I error) to test deviations
#' from Hardy–Weinberg equilibrium (HWE). Usually, \code{0.05} is used.
#' Default don't test for HWE.
#' @param only.control Use only controls (not cases) to test HWE and to
#' compute MAFs.
#' @param row.cr.min Minimun individuals' call rate that is allowed.
#' @param col.cr.min Minimum SNPs' call rate that is allowed.
#' @param maf.min Minimum Minor Allele Frequency that is allowed.
#' Usually, \code{0.01} is used. Default removes SNPs that have a zero MAF.
#' @rdname impute-qc-sub
#' @export
QC <- function(x, row.cr.min,
               col.cr.min,
               hwe.pval = NULL,
               maf.min = NULL,
               only.control = TRUE) {
  if (class(x) != "bigSNP") stop("x must be a bigSNP")

  ### HWE
  hwe.maf.qc <- function(counts) {
    if (only.control) {
      observed <- counts[4:6, ]
    } else {
      observed <- counts[1:3, ] + counts[4:6, ]
    }
    n <- colSums(observed)
    q <- (observed[1, ] + observed[2,] / 2) / n
    p <- 1 - q
    expected <- n * rbind(q^2, 2*p*q, p^2)
    bias = TRUE
    if (bias) {
      X2 <- colSums((abs(observed - expected) - 0.5)^2 / expected)
    } else {
      X2 <- colSums((observed - expected)^2 / expected)
    }
    pX2 <- stats::pchisq(X2, 1, lower.tail = F)

    maf <- pmin(p, q)

    return(list(ind.hwe.qc = which(pX2 < hwe.pval),
                ind.maf.qc = which(maf < maf.min | maf == 0)))
  }

  tmp <- Counts(x)
  tmp2 <- hwe.maf.qc(tmp$counts.col)


  ### NA COL
  n <- nrow(x$genotypes)
  call.rate.col <- colSums(tmp$counts.col) / n
  ind.cr.col.qc <- which(call.rate.col < col.cr.min)

  ### NA ROW
  m <- ncol(x$genotypes)
  call.rate.row <- 1 - tmp$counts.row / m
  ind.cr.row.qc <- which(call.rate.row < row.cr.min)

  ### Regroup everything
  ind.qc.col <- c(unlist(tmp2), ind.cr.col.qc)
  ind.qc.col <- c(ind.cr.row.qc)

  printf("Not yet implemented\n")
}

#'@description \code{sub.bigSNP}: a function
#'to get a subset of an object of class \code{bigSNP}.
#'@param ind.row Indices of the rows (individuals) to keep.
#'Negative indices can be used to exclude row indices.
#'Default: keep them all.
#'@param ind.col Indices of the columns (SNPs) to keep.
#'Negative indices can be used to exclude column indices.
#'Default: keep them all.
#'@export
#'@name sub.bigSNP
#'@rdname impute-qc-sub
sub.bigSNP <- function(x, ind.row = NULL, ind.col = NULL) {
  if (class(x) != "bigSNP") stop("x must be a bigSNP")

  if (is.null(ind.row)) ind.row <- 1:nrow(x$genotypes)
  if (is.null(ind.col)) ind.col <- 1:ncol(x$genotypes)

  newfile <- checkFile(x, "sub")
  X2 <- bigmemory::deepcopy(x$genotypes,
                            rows = ind.row,
                            cols = ind.col,
                            type = "char",
                            backingfile = paste0(newfile, ".bk"),
                            backingpath = x$backingpath,
                            descriptorfile = paste0(newfile, ".desc"))

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
#'@rdname impute-qc-sub
Impute <- function(x, ncores = 1) {
  if (class(x) != "bigSNP") stop("x must be a bigSNP")

  # get descriptors
  X.desc <- describe(x$genotypes)

  newfile <- checkFile(x, "impute")
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
