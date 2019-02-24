################################################################################

# use 2 characters for chromosomes: 1 -> 01
format_snp_id <- function(snp_id) {
  gsubfn::strapply(
    X = snp_id,
    pattern = "^(.+?)(_.+_.+_.+)$",
    FUN = function(x, y) paste0(ifelse(nchar(x) == 1, paste0("0", x), x), y),
    empty = stop("Wrong format of some SNPs."),
    simplify = 'c'
  )
}

################################################################################

#' Read variant info from one BGI file
#'
#' @param bgifile Path to one file with extension ".bgi".
#' @param snp_id Character vector of SNP IDs. These should be in the form
#'  `"<chr>_<pos>_<a1>_<a2>"` (e.g. `"1_88169_C_T"` or `"01_88169_C_T"`).
#'  **This function assumes that these IDs are uniquely identifying variants.**
#'
#' @return A data frame containing variant infos.
#'
#' @importFrom magrittr %>%
#'
#' @export
snp_readBGI <- function(bgifile, snp_id) {

  # check for packages
  if (!requireNamespace("RSQLite", quietly = TRUE))
    stop2("Please install package 'RSQLite'.")
  if (!requireNamespace("dbplyr", quietly = TRUE))
    stop2("Please install package 'dbplyr'.")

  # read variant information from index files
  db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgifile)
  on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
  infos <- dplyr::tbl(db_con, "Variant") %>%
    dplyr::mutate(
      myid = paste(chromosome, position, allele1, allele2, sep = "_")) %>%
    dplyr::filter(myid %in% snp_id) %>%
    dplyr::collect()

  # check
  ind <- match(snp_id, infos$myid)
  if (anyNA(ind)) {
    saveRDS(snp_id[is.na(ind)],
            tmp <- sub("\\.bgen\\.bgi$", "_not_found.rds", bgifile))
    stop2("Some variants have not been found (stored in '%s').", tmp)
  }

  infos
}

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
#' @param list_snp_id List (same length as the number of BGEN files) of
#'  character vector of SNP IDs to read. These should be in the form
#'  `"<chr>_<pos>_<a1>_<a2>"` (e.g. `"1_88169_C_T"` or `"01_88169_C_T"`).
#'  **This function assumes that these IDs are uniquely identifying variants.**
#' @param bgi_dir Directory of index files. Default is the same as `bgenfiles`.
#' @param ind_row An optional vector of the row indices (individuals) that
#'   are used. If not specified, all rows are used.\cr
#'   **Don't use negative indices.**
#' @param ncores Number of cores used. Default doesn't use parallelism.
#'   You may use [nb_cores()].
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
#' @examples
#' # See e.g. https://github.com/privefl/UKBiobank/blob/master/10-get-dosages.R
#'
#' @export
snp_readBGEN <- function(bgenfiles, backingfile, list_snp_id,
                         ind_row = NULL,
                         bgi_dir = dirname(bgenfiles),
                         ncores = 1) {

  # Check if backingfile already exists
  backingfile <- path.expand(backingfile)
  assert_noexist(paste0(backingfile, ".bk"))

  # Check extension of files
  sapply(bgenfiles, assert_ext, ext = "bgen")
  # Check if all files exist
  bgifiles <- file.path(bgi_dir, paste0(basename(bgenfiles), ".bgi"))
  sapply(c(bgenfiles, bgifiles), assert_exist)

  # Check list_snp_id
  assert_class(list_snp_id, "list")
  sapply(list_snp_id, assert_nona)
  assert_lengths(list_snp_id, bgenfiles)
  sizes <- lengths(list_snp_id)

  # Samples
  assert_nona(ind_row)
  N <- readBin(bgenfiles[1], what = 1L, size = 4, n = 4)[4]
  ind.row <- `if`(is.null(ind_row), seq_len(N),
                  match(seq_len(N), ind_row, nomatch = 0)) - 1L

  # Prepare Filebacked Big Matrix
  G <- FBM.code256(
    nrow = sum(ind.row >= 0),
    ncol = sum(sizes),
    code = CODE_DOSAGE,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )

  if (ncores == 1) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  # cleanup if error
  snp.info <- tryCatch(error = function(e) { unlink(G$backingfile); stop(e) }, {

    # Fill the FBM from BGEN files (and get SNP info)
    foreach(ic = seq_along(bgenfiles), .combine = 'rbind') %dopar% {

      snp_id <- format_snp_id(list_snp_id[[ic]])
      infos <- snp_readBGI(bgifiles[ic], snp_id)

      # Get dosages in FBM
      ind2 <- match(infos$myid, snp_id)
      ind.col <- (sum(sizes[seq_len(ic - 1)]) + seq_len(sizes[ic]))[ind2]
      ID <- read_bgen(
        filename = bgenfiles[ic],
        offsets = as.double(infos$file_start_position),
        BM = G,
        ind_row = ind.row,
        ind_col = ind.col,
        decode = as.raw(207 - round(0:510 * 100 / 255))
      )

      # Return variant info
      dplyr::bind_cols(infos, marker.ID = ID)[match(snp_id, infos$myid), ] %>%
        dplyr::select(chromosome, marker.ID, rsid, physical.pos = position,
                      allele1, allele2)
    }
  })

  # Create the bigSNP object
  snp.list <- structure(list(genotypes = G, map = snp.info), class = "bigSNP")

  # save it and return the path of the saved object
  rds <- sub("\\.bk$", ".rds", G$backingfile)
  saveRDS(snp.list, rds)
  rds
}

################################################################################
