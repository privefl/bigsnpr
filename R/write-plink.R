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

  # get other files
  bimfile <- sub_bed(bedfile, ".bim")
  famfile <- sub_bed(bedfile, ".fam")
  # check if files already exist
  sapply(c(bedfile, bimfile, famfile), assert_noexist)

  # create directory if doesn't exist
  assert_dir(dirname(bedfile))

  # prepare subsets
  new_fam <- x$fam[ind.row, NAMES.FAM]
  new_map <- x$map[ind.col, NAMES.MAP]
  G.round <- G$copy(code = replace(round(G$code256), is.na(G$code256), 3))
  stopifnot(all(G.round$code256 %in% 0:3))

  ## write files
  writebina(path.expand(bedfile), G.round, getInverseCode(), ind.row, ind.col)
  write.table2(new_fam, file = famfile, na = 0)
  write.table2(new_map, file = bimfile, na = 0)

  bedfile
}

################################################################################
