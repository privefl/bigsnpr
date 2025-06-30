################################################################################

#' Read PLINK files into a "bigSNP"
#'
#' Functions to read bed/bim/fam files into a [bigSNP][bigSNP-class].
#'
#' For more information on these formats, please visit
#' \href{https://www.cog-genomics.org/plink/1.9/formats#bed}{PLINK webpage}.
#' For other formats, please use PLINK to convert them in bedfiles,
#' which require minimal space to store and are faster to read. For example,
#' to convert from a VCF file, use the `--vcf` option. See [snp_plinkQC].
#'
#' @param bedfile Path to file with extension ".bed" to read.
#' You need the corresponding ".bim" and ".fam" in the same directory.
#' @param backingfile The path (without extension) for the backing files
#' for the cache of the [bigSNP][bigSNP-class] object. Default takes the bedfile
#' without the ".bed" extension.
#'
#' @return The path to the RDS file that stores the `bigSNP` object.
#' Note that this function creates another file (*.bk*) that stores the values
#' of the Filebacked Big Matrix.\cr
#' __You shouldn't read from PLINK files more than once.__ Instead, use
#' [snp_attach] to load the "bigSNP" object in any R session from backing files.
#'
#' @example examples/example-readplink.R
#' @export
snp_readBed <- function(bedfile, backingfile = sub_bed(bedfile)) {

  # Check if backingfile already exists
  backingfile <- path.expand(backingfile)
  assert_noexist(paste0(backingfile, ".bk"))

  # Get other files
  bimfile <- sub_bed(bedfile, ".bim")
  famfile <- sub_bed(bedfile, ".fam")
  # Check if all three files exist
  sapply(c(bedfile, bimfile, famfile), assert_exist)

  # Read map and family files
  fam <- bigreadr::fread2(famfile, col.names = NAMES.FAM, nThread = 1)
  bim <- bigreadr::fread2(bimfile, col.names = NAMES.MAP, nThread = 1)

  # Prepare Filebacked Big Matrix
  bigGeno <- FBM.code256(
    nrow = nrow(fam),
    ncol = nrow(bim),
    code = CODE_012,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )

  # Fill the FBM from bedfile
  reach.eof <- readbina(path.expand(bedfile), bigGeno, getCode())
  if (!reach.eof) warning("EOF of bedfile has not been reached.")

  # Create the bigSNP object
  snp.list <- structure(list(genotypes = bigGeno, fam = fam, map = bim),
                        class = "bigSNP")

  # save it and return the path of the saved object
  rds <- sub_bk(bigGeno$backingfile, ".rds")
  saveRDS(snp.list, rds)
  rds
}

################################################################################

#' @inheritParams bigsnpr-package
#' @rdname snp_readBed
#' @export
snp_readBed2 <- function(bedfile, backingfile = sub_bed(bedfile),
                         ind.row = rows_along(obj.bed),
                         ind.col = cols_along(obj.bed),
                         ncores = 1) {

  # Get mapping of bed
  obj.bed <- bed(bedfile)

  check_args()

  # Check if backingfile already exists
  backingfile <- path.expand(backingfile)
  assert_noexist(paste0(backingfile, ".bk"))

  # Read map and family files
  fam <- obj.bed$fam[ind.row, ]; rownames(fam) <- rows_along(fam)
  bim <- obj.bed$map[ind.col, ]; rownames(bim) <- rows_along(bim)

  # Prepare Filebacked Big Matrix
  bigGeno <- FBM.code256(
    nrow = nrow(fam),
    ncol = nrow(bim),
    code = CODE_012,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )

  # Fill the FBM from bedfile
  readbina2(bigGeno, obj.bed, ind.row, ind.col, ncores)

  # Create the bigSNP object
  snp.list <- structure(list(genotypes = bigGeno, fam = fam, map = bim),
                        class = "bigSNP")

  # save it and return the path of the saved object
  rds <- sub_bk(bigGeno$backingfile, ".rds")
  saveRDS(snp.list, rds)
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

  assert_exist(rdsfile)
  rdsfile <- normalizePath(rdsfile)
  assert_exist(bkfile <- sub("\\.rds$", ".bk", rdsfile))

  snp.list <- readRDS(rdsfile)
  snp.list$genotypes$backingfile <- bkfile  # in case of moving files
  snp.list$genotypes <- bigstatsr:::reconstruct_if_old(
    snp.list$genotypes, msg2 = "You should use `snp_save()`.")
  snp.list
}

################################################################################

#' Attach a "bigSNP" for examples and tests
#'
#' @param bedfile Name of one example bed file. Either
#'   - `"example.bed"` (the default),
#'   - `"example-missing.bed"`.
#'
#' @return The example "bigSNP", filebacked in the "/tmp/" directory.
#'
#' @export
snp_attachExtdata <- function(bedfile = c("example.bed", "example-missing.bed")) {

  bedfile <- system.file("extdata", match.arg(bedfile), package = "bigsnpr")
  snp_attach(
    snp_readBed(bedfile, backingfile = tempfile())
  )
}

################################################################################
