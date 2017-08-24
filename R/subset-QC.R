################################################################################

#' Subset
#'
#' A method to get a subset (copy) of a `bigSNP`.
#' This new `bigSNP` will also be backed by files (in the same directory).
#'
#' @inheritParams bigsnpr-package
#' @param ind.row Indices of the rows (individuals) to keep.
#' Negative indices __can__ be used to exclude row indices.
#' Default: keep them all.
#' @param ind.col Indices of the columns (SNPs) to keep.
#' Negative indices __can__ be used to exclude column indices.
#' Default: keep them all.
#' @param ... Not used.
#'
#' @export
#' @return The path to the RDS file that stores the `bigSNP` object.\cr
#' Note that this function creates two files whose names have been
#' automatically chosen by appending "_sub" and a number to the prefix of the
#' input bigSNP backing files.
#' @seealso [bigSNP][bigSNP-class]
#' @examples
#' str(test <- snp_attachExtdata())
#'
#' # keep only first 50 samples and SNPs
#' rdsfile <- subset(test, ind.row = 1:50, ind.col = 1:50)
#' str(snp_attach(rdsfile))
#'
#' # remove only first 50 samples and SNPs
#' rdsfile2 <- subset(test, ind.row = -(1:50), ind.col = -(1:50))
#' str(snp_attach(rdsfile2))
#'
subset.bigSNP <- function(x,
                          ind.row = rows_along(G),
                          ind.col = cols_along(G),
                          ...) {

  G <- x$genotypes
  # Support for negative indices
  ind.row <- rows_along(G)[ind.row]
  ind.col <- cols_along(G)[ind.col]

  check_args()

  # Create new FBM and fill it
  G2 <- FBM.code256(
    nrow = length(ind.row),
    ncol = length(ind.col),
    code = G$code256,
    init = NULL,
    backingfile = getNewFile(x, "sub"),
    create_bk = TRUE,
    save = FALSE
  )
  replaceSNP(G2, G, rowInd = ind.row, colInd = ind.col)

  # http://stackoverflow.com/q/19565621/6103040
  newfam <- x$fam[ind.row, , drop = FALSE]
  rownames(newfam) <- rows_along(newfam)
  newmap <- x$map[ind.col, , drop = FALSE]
  rownames(newmap) <- rows_along(newmap)

  # Create the bigSNP object
  snp.list <- structure(list(genotypes = G2,
                             fam = newfam,
                             map = newmap),
                        class = "bigSNP")

  # save it and return the path of the saved object
  rds <- sub("\\.bk$", ".rds", G2$backingfile)
  saveRDS(snp.list, rds)
  rds
}

################################################################################

# #' @name QC
# #' @description `QC`: Quality control (filters)
# #' for a `bigSNP` resulting in a `bigSNP` of lower dimension.
# #' @param hwe.pval Level threshold (allowed type-I error) to test deviations
# #' from Hardy–Weinberg equilibrium (HWE) from controls only. Default is `1e-6`.
# #' @param rm.sex Keep only SNPs on the first 22 chrosmosomes? Default is `FALSE`.
# #' @param row.cr.min Minimun individuals' call rate that is allowed.
# #' Default is 95\%.
# #' @param col.cr.min Minimum SNPs' call rate that is allowed. Default is 95\%.
# #' @param maf.min Minimum Minor Allele Frequency that is allowed.
# #' Usually, `0.01` is used. Default only removes SNPs that have a zero MAF.
# #' @rdname impute-qc-sub
# #' @export
# QC <- function(x, row.cr.min = 0.95,
#                col.cr.min = 0.95,
#                hwe.pval = 1e-6,
#                maf.min = NULL,
#                rm.sex = FALSE) {
#   check_x(x)
#
#   counts <- snp_counts(x)
#
#   ### HWE
#   hwe.qc <- function(observed) {
#     n <- colSums(observed)
#     q <- (observed[1, ] + observed[2,] / 2) / n
#     p <- 1 - q
#     expected <- n * rbind(q^2, 2*p*q, p^2)
#
#     #X2 <- colSums((abs(observed - expected) - 0.5)^2 / expected)
#     X2 <- colSums((observed - expected)^2 / expected)
#     pX2 <- stats::pchisq(X2, 1, lower.tail = FALSE)
#
#     return(which(pX2 < hwe.pval))
#   }
#   ind.hwe.qc <- hwe.qc(counts$cols.controls) # only controls
#
#   ### MAF
#   # controls + cases
#   observed <- counts$cols.controls + counts$cols.cases
#   n <- colSums(observed)
#   q <- (observed[1, ] + observed[2,] / 2) / n
#   maf <- pmin(q, 1 - q)
#   ind.maf.qc <- which(maf < maf.min | maf == 0)
#
#   ### NA COL
#   n.all <- nrow(x$genotypes)
#   call.rate.col <- n / n.all
#   ind.cr.col.qc <- which(call.rate.col < col.cr.min)
#
#   ### NOT AUTOSOMAL
#   if (rm.sex) {
#     ind.sex <- which(x$map$chromosome > 22)
#   } else {
#     ind.sex <- integer(0)
#   }
#
#   ### NA ROW
#   m.all <- ncol(x$genotypes)
#   call.rate.row <- 1 - counts$rows / m.all
#   ind.cr.row.qc <- which(call.rate.row < row.cr.min)
#
#   ### Regroup everything
#   ind.qc.col <- c(ind.hwe.qc, ind.maf.qc, ind.cr.col.qc, ind.sex)
#   ind.qc.row <- c(ind.cr.row.qc)
#
#   return(sub.bigSNP(x,
#                     ind.row = `if`(length(ind.qc.row) > 0,
#                                    -ind.qc.row, seq(n.all)),
#                     ind.col = `if`(length(ind.qc.col) > 0,
#                                    -ind.qc.col, seq(m.all))))
# }

################################################################################

# fisher2 <- function(nbNA.ca, nbNA, N.ca, N, relErr = 1 + 1e-7) {
#   # precomputation
#   dhyper.precomputed <- vector("list", max(nbNA))
#   for (m in sort(unique(nbNA))) {
#     n <- N - m
#     lo <- max(0L, N.ca - n)
#     hi <- min(N.ca, m)
#     d <- dhyper(lo:hi, m, n, N.ca)
#     dhyper.precomputed[[m]] <- d / sum(d)
#   }
#
#   # compute p-values
#   len <- length(nbNA)
#   PVAL <- numeric(len)
#   ind <- nbNA.ca - pmax(0L, N.ca - N + nbNA) + 1
#   for (i in seq_len(len)) {
#     d <- dhyper.precomputed[[nbNA[i]]]
#     PVAL[i] <- sum(d[d <= (d[ind[i]] * relErr)])
#   }
#
#   PVAL
# }

################################################################################
