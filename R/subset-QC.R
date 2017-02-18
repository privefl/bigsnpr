################################################################################

#' @title Imputation, quality control and subset for a "bigSNP".
#' @name impute-qc-sub
#' @inheritParams bigsnpr-package
#' @return A new `bigSNP`.
#' @seealso [bigSNP][bigSNP-class]
#' @example examples/example.sub.bigSNP.R
NULL


################################################################################

#' Copy (subset) of a "bigSNP"
#'
#' @description `sub.bigSNP`: a function
#' to get a subset of an object of class `bigSNP`.
#' @param ind.row Indices of the rows (individuals) to keep.
#' Negative indices can be used to exclude row indices.
#' Default: keep them all.
#' @param ind.col Indices of the columns (SNPs) to keep.
#' Negative indices can be used to exclude column indices.
#' Default: keep them all.
#' @param backed Should the new `bigSNP` be filebacked? Default is `TRUE`.
#' @param shared Should the new genotype matrix be shared? Default is `TRUE`.
#' @export
#' @name sub.bigSNP
#' @rdname impute-qc-sub
sub.bigSNP <- function(x, ind.row = seq(nrow(x$genotypes)),
                       ind.col = seq(ncol(x$genotypes)),
                       backed = TRUE, shared = TRUE) {
  check_x(x)

  if (backed) {
    newfile <- checkFile(x, "sub")
    X2 <- deepcopy(x$genotypes,
                   rows = ind.row,
                   cols = ind.col,
                   # type = "char",
                   backingfile = paste0(newfile, ".bk"),
                   backingpath = x$backingpath,
                   descriptorfile = paste0(newfile, ".desc"))

    # http://stackoverflow.com/q/19565621/6103040
    newfam <- x$fam[ind.row, ]
    rownames(newfam) <- seq_len(nrow(newfam))
    newmap <- x$map[ind.col, ]
    rownames(newmap) <- seq_len(nrow(newmap))

    snp_list <- list(genotypes = X2,
                     fam = newfam,
                     map = newmap,
                     backingfile = newfile,
                     backingpath = x$backingpath)
    class(snp_list) <- "bigSNP"

    saveRDS(snp_list, file.path(x$backingpath, paste0(newfile, ".rds")))

  } else {
    X2 <- deepcopy(x$genotypes,
                   rows = ind.row,
                   cols = ind.col,
                   # type = "char",
                   shared = shared)

    snp_list <- list(genotypes = X2,
                     fam = x$fam[ind.row, ],
                     map = x$map[ind.col, ],
                     backingfile = NULL,
                     backingpath = NULL)
    class(snp_list) <- "bigSNP"
  }

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

  counts <- snp_counts(x)

  ### HWE
  hwe.qc <- function(observed) {
    n <- colSums(observed)
    q <- (observed[1, ] + observed[2,] / 2) / n
    p <- 1 - q
    expected <- n * rbind(q^2, 2*p*q, p^2)

    #X2 <- colSums((abs(observed - expected) - 0.5)^2 / expected)
    X2 <- colSums((observed - expected)^2 / expected)
    pX2 <- stats::pchisq(X2, 1, lower.tail = FALSE)

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
