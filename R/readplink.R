################################################################################

#' Read PLINK files into a "bigSNP"
#'
#' Function to read bed/bim/fam files into a [bigSNP][bigSNP-class].
#'
#' For more information on these formats, please visit
#' \href{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed}{PLINK webpage}.
#' For other formats, please use PLINK to convert them in bedfiles,
#' which require minimal space to store and are faster to read.
#' For example, to convert from a VCF file, use the `--vcf` option.
#'
#' @param bedfile Path to file with extension ".bed" to read.
#' You need the corresponding ".bim" and ".fam" in the same directory.
#' @param backingfile The root name for the backing file(s) for the cache of
#' the [bigSNP][bigSNP-class] object.
#' @param backingpath The path to the directory containing the backing files.
#' Default is `"backingfiles"`.
#'
#' @return The path to the RDS file that stores the `bigSNP` object.
#' Note that this function creates two other files which store the
#' `big.matrix` and its descriptor.\cr
#' __You shouldn't read from PLINK files more than once.__ Instead, use
#' [snp_attach] to load the "bigSNP" object in any R session from backing files.
#'
#' @example examples/example-readplink.R
#' @export
snp_readBed <- function(bedfile, backingfile,
                        backingpath = "backingfiles") {

  # create backingpath if doesn't exist
  backingpath <- path.expand(backingpath)
  assert_dir(backingpath)

  # check if backingfile already exists
  rootPath <- file.path(backingpath, backingfile)
  assert_noexist(paste0(rootPath, ".bk"))

  # check extension of file
  assert_ext(bedfile, "bed")
  # get other files
  bimfile <- sub("\\.bed$", ".bim", bedfile)
  famfile <- sub("\\.bed$", ".fam", bedfile)
  # check if all three files exist
  sapply(c(bedfile, bimfile, famfile), assert_exist)

  # read map and family files
  fam <- data.table::fread(famfile, data.table = FALSE)
  names(fam) <- NAMES.FAM
  bim <- data.table::fread(bimfile, data.table = FALSE)
  names(bim) <- NAMES.MAP

  # prepare big.matrix
  n <- nrow(fam)
  m <- nrow(bim)
  bigGeno <- big.matrix(n, m, type = "raw",
                        backingfile = paste0(backingfile, ".bk"),
                        backingpath = backingpath,
                        descriptorfile = paste0(backingfile, ".desc"))
  # removing unnecessary ".desc" file
  file.remove(paste0(rootPath, ".desc"))

  # fill the `big.matrix` from bedfile
  reach.eof <- readbina(bedfile, bigGeno@address, getCode())
  if (!reach.eof) warning("EOF of bedfile has not been reached.")

  bigGeno.code <- as.BM.code(bigGeno, code = CODE_012)

  # create the `bigSNP`
  rds <- paste0(rootPath, ".rds")
  snp_list <- structure(list(genotypes = describe(bigGeno.code),
                             fam = fam,
                             map = bim,
                             savedIn = rds),
                        class = "bigSNP")

  # save it and return the path of the saved object
  saveRDS(snp_list, rds)
  rds
}

################################################################################

#' Attach a "bigSNP" from backing files
#'
#' Load a [bigSNP][bigSNP-class] from backing files into R.
#'
#' This is often just a call to [readRDS]. But it also checks if you have moved
#' the two (".bk" and ".rds") backing files to another directory.
#'
#' @param rdsfile The path of the ".rds" which stores the `bigSNP` object.
#'
#' @return The `bigSNP` object.
#' @example examples/example-readplink.R
#'
#' @export
snp_attach <- function(rdsfile) {

  rdsfile <- path.expand(rdsfile)
  snp.list <- readRDS(rdsfile)

  if (rdsfile != snp.list$savedIn) {
    # check if backing file exists
    bkfile <- sub("\\.rds$", ".bk", rdsfile)
    assert_exist(bkfile)
    # change relative paths
    snp.list$genotypes@description$dirname <- dirname(bkfile)
    snp.list$savedIn <- rdsfile
  }

  snp.list
}

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
                         ind.row = rows_along(X),
                         ind.col = cols_along(X)) {

  X <- attach.BM(x$genotypes)

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
  X@code <- replace(round(X@code), is.na(X@code), 3)
  stopifnot(all(X@code %in% 0:3))

  writebina(bedfile, X, getInverseCode(), ind.row, ind.col)

  bedfile
}

################################################################################

#' Attach a "bigSNP" for examples and tests
#'
#' @return The example "bigSNP", filebacked in the "/tmp/" directory.
#'
#' @export
snp_attachExtdata <- function() {

  tmpfile <- tempfile()
  bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
  rdsfile <- snp_readBed(bedfile,
                         backingfile = basename(tmpfile),
                         backingpath = dirname(tmpfile))

  snp_attach(rdsfile)
}

################################################################################
