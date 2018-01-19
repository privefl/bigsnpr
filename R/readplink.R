################################################################################

#' Read PLINK files into a "bigSNP"
#'
#' Function to read bed/bim/fam files into a [bigSNP][bigSNP-class].
#'
#' For more information on these formats, please visit
#' \href{https://www.cog-genomics.org/plink/1.9/formats#bed}{PLINK webpage}.
#' For other formats, please use PLINK to convert them in bedfiles,
#' which require minimal space to store and are faster to read. For example,
#' to convert from a VCF file, use the `--vcf` option. See [snp_plinkQC].
#'
#' @param bedfile Path to file with extension ".bed" to read.
#' You need the corresponding ".bim" and ".fam" in the same directory.
#' @param backingfile The path (without extension) for the backing file(s)
#' for the cache of the [bigSNP][bigSNP-class] object. Default takes the bedfile
#' without the ".bed" extension.
#'
#' @return The path to the RDS file that stores the `bigSNP` object.
#' Note that this function creates one other file which stores the values of
#' the Filebacked Big Matrix.\cr
#' __You shouldn't read from PLINK files more than once.__ Instead, use
#' [snp_attach] to load the "bigSNP" object in any R session from backing files.
#'
#' @example examples/example-readplink.R
#' @export
snp_readBed <- function(bedfile, backingfile = sub("\\.bed$", "", bedfile)) {

  # Check if backingfile already exists
  backingfile <- path.expand(backingfile)
  assert_noexist(paste0(backingfile, ".bk"))

  # Check extension of file
  assert_ext(bedfile, "bed")
  # Get other files
  bimfile <- sub("\\.bed$", ".bim", bedfile)
  famfile <- sub("\\.bed$", ".fam", bedfile)
  # Check if all three files exist
  sapply(c(bedfile, bimfile, famfile), assert_exist)

  # Read map and family files
  fam <- data.table::fread(famfile, data.table = FALSE)
  names(fam) <- NAMES.FAM
  bim <- data.table::fread(bimfile, data.table = FALSE)
  names(bim) <- NAMES.MAP

  # Prepare Filebacked Big Matrix
  bigGeno <- FBM.code256(
    nrow = nrow(fam),
    ncol = nrow(bim),
    code = CODE_012,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE,
    save = FALSE
  )

  # Fill the FBM from bedfile
  reach.eof <- readbina(bedfile, bigGeno, getCode())
  if (!reach.eof) warning("EOF of bedfile has not been reached.")

  # Create the bigSNP object
  snp.list <- structure(list(genotypes = bigGeno,
                             fam = fam,
                             map = bim),
                        class = "bigSNP")

  # save it and return the path of the saved object
  rds <- sub("\\.bk$", ".rds", bigGeno$backingfile)
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

  rdsfile <- normalizePath(rdsfile)
  assert_exist(bkfile <- sub("\\.rds$", ".bk", rdsfile))

  snp.list <- readRDS(rdsfile)
  snp.list$genotypes$backingfile <- bkfile  # in case of moving files
  snp.list
}

################################################################################

#' Attach a "bigSNP" for examples and tests
#'
#' @inheritParams snp_readBed
#'
#' @return The example "bigSNP", filebacked in the "/tmp/" directory.
#'
#' @export
snp_attachExtdata <- function(bedfile = system.file("extdata", "example.bed",
                                                    package = "bigsnpr")) {

  snp_attach(
    snp_readBed(bedfile, backingfile = tempfile())
  )
}

################################################################################
