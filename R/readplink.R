################################################################################

#' Replace extension '.bed'
#'
#' @param path String with extension '.bed'.
#' @param replacement Replacement of '.bed'. Default replaces by nothing.
#'   Can be useful to replace e.g. by '.bim' or '.fam'.
#' @param stop_if_not_ext If `replacement != ""`, whether to error if
#'   replacement is not an extension (starting with a '.').
#'
#' @return String with extension '.bed' replaced by `replacement`.
#' @export
#'
#' @examples
#' path <- "toto.bed"
#' sub_bed(path)
#' sub_bed(path, ".bim")
#' sub_bed(path, ".fam")
#' sub_bed(path, "_QC", stop_if_not_ext = FALSE)
sub_bed <- function(path, replacement = "", stop_if_not_ext = TRUE) {
  pattern <- "\\.bed$"
  if (!grepl(pattern, path))
    stop2("Path '%s' must have 'bed' extension.", path)
  if (stop_if_not_ext && nchar(replacement) > 0 && substr(replacement, 1, 1) != ".")
    stop2("Replacement must be an extension starting with '.' if provided.")
  sub(pattern, replacement, path)
}

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
#' @param backingfile The path (without extension) for the backing files
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
  fam <- bigreadr::fread2(famfile, col.names = NAMES.FAM)
  bim <- bigreadr::fread2(bimfile, col.names = NAMES.MAP)

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
  reach.eof <- readbina(bedfile, bigGeno, getCode())
  if (!reach.eof) warning("EOF of bedfile has not been reached.")

  # Create the bigSNP object
  snp.list <- structure(list(genotypes = bigGeno,
                             fam = fam,
                             map = bim),
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

  rdsfile <- normalizePath(rdsfile)
  assert_exist(bkfile <- sub("\\.rds$", ".bk", rdsfile))

  snp.list <- readRDS(rdsfile)
  snp.list$genotypes$backingfile <- bkfile  # in case of moving files
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
snp_attachExtdata <- function(bedfile = "example.bed") {

  bedfile <- system.file("extdata", bedfile, package = "bigsnpr")
  snp_attach(
    snp_readBed(bedfile, backingfile = tempfile())
  )
}

################################################################################
