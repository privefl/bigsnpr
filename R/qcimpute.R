################################################################################

#' @title Imputation, quality control and subset for a "bigSNP".
#' @name impute-qc-sub
#' @inheritParams bigsnpr-package
#' @return A new `bigSNP`.
#' @seealso [bigSNP][bigSNP-class]
#' @example examples/example.sub.bigSNP.R
NULL


################################################################################

#' @description `sub.bigSNP`: a function
#' to get a subset of an object of class `bigSNP`.
#' @param ind.row Indices of the rows (individuals) to keep.
#' Negative indices can be used to exclude row indices.
#' Default: keep them all.
#' @param ind.col Indices of the columns (SNPs) to keep.
#' Negative indices can be used to exclude column indices.
#' Default: keep them all.
#' @export
#' @name sub.bigSNP
#' @rdname impute-qc-sub
sub.bigSNP <- function(x, ind.row = seq(nrow(x$genotypes)),
                       ind.col = seq(ncol(x$genotypes))) {
  check_x(x)

  newfile <- checkFile(x, "sub")
  X2 <- deepcopy(x$genotypes,
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

#' @name QC
#' @description `QC`: Quality control (filters)
#' for a `bigSNP` resulting in a `bigSNP` of lower dimension.
#' @param hwe.pval Level threshold (allowed type-I error) to test deviations
#' from Hardy–Weinberg equilibrium (HWE) from controls only. Default is `1e-6`.
#' @param rm.sex Keep only SNPs on the first 22 chrosmosomes? Default is `FALSE`.
#' @param row.cr.min Minimun individuals' call rate that is allowed.
#' Default is 95\%.
#' @param col.cr.min Minimum SNPs' call rate that is allowed. Default is 95\%.
#' @param maf.min Minimum Minor Allele Frequency that is allowed.
#' Usually, `0.01` is used. Default only removes SNPs that have a zero MAF.
#' @rdname impute-qc-sub
#' @export
QC <- function(x, row.cr.min = 0.95,
               col.cr.min = 0.95,
               hwe.pval = 1e-6,
               maf.min = NULL,
               rm.sex = FALSE) {
  check_x(x)

  counts <- Counts(x)

  ### HWE
  hwe.qc <- function(observed) {
    n <- colSums(observed)
    q <- (observed[1, ] + observed[2,] / 2) / n
    p <- 1 - q
    expected <- n * rbind(q^2, 2*p*q, p^2)

    #X2 <- colSums((abs(observed - expected) - 0.5)^2 / expected)
    X2 <- colSums((observed - expected)^2 / expected)
    pX2 <- stats::pchisq(X2, 1, lower.tail = F)

    return(which(pX2 < hwe.pval))
  }
  ind.hwe.qc <- hwe.qc(counts$cols.controls) # only controls

  ### MAF
  # controls + cases
  observed <- counts$cols.controls + counts$cols.cases
  n <- colSums(observed)
  q <- (observed[1, ] + observed[2,] / 2) / n
  maf <- pmin(q, 1 - q)
  ind.maf.qc <- which(maf < maf.min | maf == 0)

  ### NA COL
  n.all <- nrow(x$genotypes)
  call.rate.col <- n / n.all
  ind.cr.col.qc <- which(call.rate.col < col.cr.min)

  ### NOT AUTOSOMAL
  if (rm.sex) {
    ind.sex <- which(x$map$chromosome > 22)
  } else {
    ind.sex <- integer(0)
  }

  ### NA ROW
  m.all <- ncol(x$genotypes)
  call.rate.row <- 1 - counts$rows / m.all
  ind.cr.row.qc <- which(call.rate.row < row.cr.min)

  ### Regroup everything
  ind.qc.col <- c(ind.hwe.qc, ind.maf.qc, ind.cr.col.qc, ind.sex)
  ind.qc.row <- c(ind.cr.row.qc)

  return(sub.bigSNP(x,
                    ind.row = `if`(length(ind.qc.row) > 0,
                                   -ind.qc.row, seq(n.all)),
                    ind.col = `if`(length(ind.qc.col) > 0,
                                   -ind.qc.col, seq(m.all))))
}

################################################################################

#'@description `Impute`: Imputation function
#'for a `bigSNP`.
#'@param verbose Print progress? Default is `FALSE`.
#'@export
#'@name Impute
#'@rdname impute-qc-sub
Impute <- function(x, ncores = 1, verbose = FALSE) {
  check_x(x)

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
    if (verbose)
      printf("Imputing chromosome %d with \"nearest neighbors\"...\n", lims[3])

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

    if (verbose)
      printf("Done imputing chromosome %d.\n", lims[3])

    return(0)
  }

  range.chr <- LimsChr(x)

  obj <- foreach::foreach(i = 1:nrow(range.chr),
                          .noexport = c("x", "X2"),
                          .packages = "bigmemory")
  expr_fun <- function(i) {
    ImputeChr(range.chr[i, ])
  }
  foreach2(obj, expr_fun, ncores,
           outfile = `if`(verbose, "", NULL))

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
