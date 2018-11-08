################################################################################

DECODE_BGEN <- as.raw(207 - round(0:510 * 100 / 255))

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
#'   The corresponding ".bgen.bgi" index files must exist.
#' @param backingfile The path (without extension) for the backing files
#'   for the cache of the [bigSNP][bigSNP-class] object.
#' @param list_snp_pos List (same length as the number of BGEN files) of
#'  SNP position to read.
#' @param bgi_dir Directory of index files. Default is the same as `bgenfiles`.
#'
#' @return The path to the RDS file that stores the `bigSNP` object.
#' Note that this function creates one other file which stores the values of
#' the Filebacked Big Matrix.\cr
#' __You shouldn't read from BGEN files more than once.__ Instead, use
#' [snp_attach] to load the "bigSNP" object in any R session from backing files.
#'
#' @importFrom magrittr %>%
#' @import foreach
#'
#' @export
snp_readBGEN <- function(bgenfiles, backingfile, list_snp_pos,
                         bgi_dir = dirname(bgenfiles)) {

  if (!requireNamespace("RSQLite", quietly = TRUE))
    stop2("Please install package 'RSQLite'.")

  # Check if backingfile already exists
  backingfile <- path.expand(backingfile)
  assert_noexist(paste0(backingfile, ".bk"))

  # Check extension of files
  sapply(bgenfiles, assert_ext, ext = "bgen")
  # Check if all files exist
  bgifiles <- file.path(bgi_dir, paste0(basename(bgenfiles), ".bgi"))
  sapply(c(bgenfiles, bgifiles), assert_exist)

  # Check list_snp_pos
  assert_class(list_snp_pos, "list")
  sapply(list_snp_pos, assert_nona)
  sizes <- lengths(list_snp_pos)
  assert_lengths(sizes, bgenfiles)

  # Prepare Filebacked Big Matrix
  G <- FBM.code256(
    nrow = readBin(bgenfiles[1], what = 1L, size = 4, n = 4)[4],
    ncol = sum(sizes),
    code = CODE_DOSAGE,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )

  # Fill the FBM from BGEN files (and get SNP info)
  snp.info <- foreach(ic = seq_along(bgenfiles), .combine = 'rbind') %do% {

    snp_pos <- list_snp_pos[[ic]]

    # Read variant info (+ position in file) from index files
    db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgifiles[ic])
    infos <- dplyr::tbl(db_con, "Variant") %>%
      dplyr::filter(position %in% snp_pos) %>%
      dplyr::collect()
    RSQLite::dbDisconnect(db_con)
    offsets <- as.double(infos$file_start_position)

    # Check if found all SNPs
    ind <- match(snp_pos, infos$position)
    if (anyNA(ind)) stop2("Some variants have not been found.")

    # Get dosages in FBM
    ind.col <- sum(sizes[seq_len(ic - 1)]) + seq_len(sizes[ic])
    read_bgen(bgenfiles, offsets, G, ind.col[match(infos$position, snp_pos)],
              DECODE_BGEN)

    # Return variant info
    stats::setNames(infos[ind, c(1, 3, 2, 5, 6)], NAMES.MAP[-3])
  }

  # Create the bigSNP object
  snp.list <- structure(
    list(genotypes = G, map = snp.info),
    class = "bigSNP"
  )

  # save it and return the path of the saved object
  rds <- sub("\\.bk$", ".rds", G$backingfile)
  saveRDS(snp.list, rds)
  rds
}

################################################################################
