################################################################################

#' Write PLINK files from a "bigSNP"
#'
#' Function to write bed/bim/fam files from a [bigSNP][bigSNP-class].
#' This will use the slot `code` **rounded** to write 0s, 1s, 2s or NAs.
#'
#' @inheritParams bigsnpr-package
#' @param bedfile Path to file with extension ".bed" to create.
#'
#' @return The input `bedfile` path.
#'
#' @example examples/example-writeplink.R
#' @export
snp_writeBed <- function(x, bedfile,
                         ind.row = rows_along(G),
                         ind.col = cols_along(G)) {

  G <- x$genotypes

  check_args(G = "assert_class(G, 'FBM.code256')")

  # check extension of file
  assert_ext(bedfile, "bed")
  # get other files
  bimfile <- sub("\\.bed$", ".bim", bedfile)
  famfile <- sub("\\.bed$", ".fam", bedfile)
  # check if files already exist
  sapply(c(bedfile, bimfile, famfile), assert_noexist)

  # create directory if doesn't exist
  assert_dir(dirname(bedfile))

  # write map and family files
  write.table2(x$fam[ind.row, NAMES.FAM], file = famfile)
  write.table2(x$map[ind.col, NAMES.MAP], file = bimfile)

  ## write bed file
  G.round <- G$copy(code = replace(round(G$code256), is.na(G$code256), 3))
  stopifnot(all(G.round$code256 %in% 0:3))

  writebina(bedfile, G.round, getInverseCode(), ind.row, ind.col)

  bedfile
}

################################################################################
