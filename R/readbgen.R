################################################################################

#' Read BGEN files into a "bigSNP"
#'
#' Function to read the UK Biobank BGEN files into a [bigSNP][bigSNP-class].
#'
#' For more information on this format, please visit
#' \href{https://bitbucket.org/gavinband/bgen/}{BGEN webpage}.
#'
#' This function is designed to read UK Biobank imputation files. This assumes
#' that variants have been compressed with zlib, that there are only 2 possible
#' alleles, and that each probability is stored on 8 bits.
#'
#' @param bgenfiles Character vector of paths to files with extension ".bgen".
#' @param backingfile The path (without extension) for the backing files
#'   for the cache of the [bigSNP][bigSNP-class] object.
#' @param snp_id Character vector of SNP IDs (not rsid) to read
#'   (e.g. "1:88169_C_T" in the UK Biobank).
#'
#' @return The path to the RDS file that stores the `bigSNP` object.
#' Note that this function creates one other file which stores the values of
#' the Filebacked Big Matrix.\cr
#' __You shouldn't read from BGEN files more than once.__ Instead, use
#' [snp_attach] to load the "bigSNP" object in any R session from backing files.
#'
#' @export
snp_readBGEN <- function(bgenfiles, backingfile, snp_id) {

  # Check if backingfile already exists
  backingfile <- path.expand(backingfile)
  assert_noexist(paste0(backingfile, ".bk"))

  # Check extension of files
  sapply(bgenfiles, assert_ext, ext = "bgen")
  # Check if all files exist
  sapply(bgenfiles, assert_exist)

  # Prepare Filebacked Big Matrix
  bigGeno <- FBM.code256(
    nrow = readBin(bgenfiles[[1]], what = 1L, size = 4, n = 4)[4],
    ncol = length(snp_id),
    code = CODE_DOSAGE,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )

  # Fill the FBM from BGEN files (and get SNP info)
  snp.info <- readbgen(bgenfiles, snp_id, bigGeno)

  # Create the bigSNP object
  snp.list <- structure(
    list(genotypes = bigGeno, map = snp.info),
    class = "bigSNP"
  )

  # save it and return the path of the saved object
  rds <- sub("\\.bk$", ".rds", bigGeno$backingfile)
  saveRDS(snp.list, rds)
  rds
}

################################################################################
