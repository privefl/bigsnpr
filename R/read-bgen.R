################################################################################

# use 2 characters for chromosomes: 1 -> 01
format_snp_id <- function(snp_id) {
  # if chr is only one digit, append a '0' to the beginning
  snp_id <- ifelse(substr(snp_id, 2, 2) == "_", paste0("0", snp_id), snp_id)
  if (any(substr(snp_id, 3, 3) != "_")) stop("Wrong format of some variants.")
  snp_id
}

################################################################################

#' Read variant info from one BGI file
#'
#' @param bgifile Path to one file with extension ".bgi".
#' @param snp_id Character vector of SNP IDs. These should be in the form
#'  `"<chr>_<pos>_<a1>_<a2>"` (e.g. `"1_88169_C_T"` or `"01_88169_C_T"`).
#'  **This function assumes that these IDs are uniquely identifying variants.**
#'  Default is `NULL`, and returns information on all variants.
#'
#' @return A data frame containing variant information.
#'
#' @importFrom magrittr %>%
#'
#' @export
snp_readBGI <- function(bgifile, snp_id = NULL) {

  # check for packages
  assert_package("RSQLite")
  assert_package("dbplyr")

  # read variant information from index files
  db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgifile)
  on.exit(RSQLite::dbDisconnect(db_con), add = TRUE, after = FALSE)

  if (is.null(snp_id)) {

    dplyr::tbl(db_con, "Variant") %>%
      dplyr::mutate(file_start_position = as.double(file_start_position)) %>%
      dplyr::collect()

  } else {

    snp_id <- format_snp_id(snp_id)
    snp_pos <- as.integer(sub("^[[:alnum:]]{2}_([[:digit:]]+)_.+$", "\\1", snp_id))
    info <- dplyr::tbl(db_con, "Variant") %>%
      dplyr::filter(position %in% snp_pos) %>%
      dplyr::mutate(file_start_position = as.double(file_start_position)) %>%
      dplyr::collect()

    # check
    info_id <- with(info, paste(chromosome, position, allele1, allele2, sep = "_"))
    ind <- match(snp_id, format_snp_id(info_id))
    if (anyNA(ind)) {
      saveRDS(snp_id[is.na(ind)],
              tmp <- sub("\\.bgen\\.bgi$", "_not_found.rds", bgifile))
      stop2("Some variants have not been found (stored in '%s').", tmp)
    }

    info[ind, ]

  }
}

################################################################################

check_bgen_format <- function(bgenfile) {

  header <- readBin(bgenfile, what = integer(), size = 4, n = 5)
  bgen_int <- readBin(charToRaw("bgen"), what = integer(), size = 4)
  if (!identical(header[5], bgen_int))
    stop2("'%s' is not a BGEN file.", bgenfile)

  header_raw <- readBin(bgenfile, what = raw(), n = 4 + header[2])
  flags <- rawToBits(tail(header_raw, 4))
  if (!identical(flags[1:2], rawToBits(as.raw(1))[1:2]))
    stop2("'%s' is not compressed with zlib.", bgenfile)
  if (!identical(flags[3:6], rawToBits(as.raw(2))[1:4]))
    stop2("'%s' is not using Layout 2.", bgenfile)

  N <- header[4]
}

################################################################################

#' Read BGEN files into a "bigSNP"
#'
#' Function to read the UK Biobank BGEN files into a [bigSNP][bigSNP-class].
#'
#' For more information on this format, please visit
#' \href{https://code.enkre.net/bgen}{BGEN webpage}.
#'
#' This function is designed to read UK Biobank imputation files. This assumes
#' that variants have been compressed with zlib, that there are only 2 possible
#' alleles, and that each probability is stored on 8 bits. For example, if you
#' use *qctool* to generate your own BGEN files, please make sure you are using
#' options '`-ofiletype bgen_v1.2 -bgen-bits 8 -assume-chromosome`'.
#'
#' If the format is not the expected one, this will result in an error or even
#' a crash of your R session. Another common source of error is due to corrupted
#' files; e.g. if using UK Biobank files, compare the result of [tools::md5sum()]
#' with the ones at \url{https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=998}.
#'
#' You can look at some example code from my papers on how to use this function:
#' - \url{https://github.com/privefl/paper-infer/blob/main/code/prepare-geno-simu.R}
#' - \url{https://github.com/privefl/paper-misspec/blob/main/code/prepare-genotypes.R}
#'
#' @param bgenfiles Character vector of paths to files with extension ".bgen".
#'   The corresponding ".bgen.bgi" index files must exist.
#' @param backingfile The path (without extension) for the backing files (".bk"
#'   and ".rds") that are created by this function for storing the
#'   [bigSNP][bigSNP-class] object.
#' @param list_snp_id List of character vectors of SNP IDs to read, with one
#'   vector per BGEN file. Each SNP ID should be in the form
#'   `"<chr>_<pos>_<a1>_<a2>"` (e.g. `"1_88169_C_T"` or `"01_88169_C_T"`).
#'   If you have one BGEN file only, just wrap your vector of IDs with `list()`.
#'   **This function assumes that these IDs are uniquely identifying variants.**
#' @param bgi_dir Directory of index files. Default is the same as `bgenfiles`.
#' @param ind_row An optional vector of the row indices (individuals) that
#'   are used. If not specified, all rows are used. **Don't use negative indices.**
#'   You can access the sample IDs corresponding to the genotypes from the *.sample*
#'   file, and use e.g. `match()` to get indices corresponding to the ones you want.
#'
#' @param ncores Number of cores used. Default doesn't use parallelism.
#'   You may use [bigstatsr::nb_cores()].
#' @param read_as How to read BGEN probabilities? Currently implemented:
#'   - as dosages (rounded to two decimal places), the default,
#'   - as hard calls, randomly sampled based on those probabilities
#'   (similar to PLINK option '`--hard-call-threshold random`').
#'
#' @return The path to the RDS file `<backingfile>.rds` that stores the `bigSNP`
#'   object created by this function.\cr
#'   Note that this function creates another file (*.bk*) that stores the values
#'   of the FBM (`$genotypes`). The rows corresponds to the order of `ind_row`;
#'   the columns to the order of `list_snp_id`. The `$map` component of the
#'   `bigSNP` object stores some information on the variants (including allele
#'   frequencies and INFO scores computed from the imputation probabilities).
#'   However, it does not have a `$fam` component; you should use the individual
#'   IDs in the *.sample* file (filtered with `ind_row`) to add external
#'   information on the individuals.\cr
#'   __You shouldn't read from BGEN files more than once.__ Instead, use
#'   [snp_attach] to load the "bigSNP" object in any R session from backing files.
#'
#' @importFrom magrittr %>%
#'
#' @export
snp_readBGEN <- function(bgenfiles, backingfile, list_snp_id,
                         ind_row = NULL,
                         bgi_dir = dirname(bgenfiles),
                         read_as = c("dosage", "random"),
                         ncores = 1) {

  dosage <- identical(match.arg(read_as), "dosage")

  # Check if backingfile already exists
  backingfile <- path.expand(backingfile)
  assert_noexist(paste0(backingfile, ".bk"))

  # Check extension of files
  bgenfiles <- path.expand(bgenfiles)
  sapply(bgenfiles, assert_ext, ext = "bgen")
  # Check if all files exist
  bgifiles <- file.path(bgi_dir, paste0(basename(bgenfiles), ".bgi"))
  sapply(c(bgenfiles, bgifiles), assert_exist)

  # Check list_snp_id
  assert_class(list_snp_id, "list")
  sapply(list_snp_id, assert_nona)
  assert_lengths(list_snp_id, bgenfiles)
  sizes <- lengths(list_snp_id)

  # Check format of BGEN files + check samples
  all_N <- sapply(bgenfiles, check_bgen_format)
  N <- all_N[1]
  bigassertr::assert_all(all_N, N)
  if (is.null(ind_row)) ind_row <- seq_len(N)
  assert_nona(ind_row)
  stopifnot(all(ind_row >= 1 & ind_row <= N))

  # Prepare Filebacked Big Matrix
  G <- FBM.code256(
    nrow = length(ind_row),
    ncol = sum(sizes),
    code = CODE_DOSAGE,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )

  # cleanup if error
  snp.info <- tryCatch(error = function(e) { unlink(G$backingfile); stop(e) }, {

    # Fill the FBM from BGEN files (and get SNP info)
    do.call("rbind", lapply(seq_along(bgenfiles), function(ic) {

      info <- snp_readBGI(bgifiles[ic], list_snp_id[[ic]])

      # Get dosages in FBM
      ind.col <- sum(sizes[seq_len(ic - 1)]) + seq_len(sizes[ic])
      varinfo <- read_bgen(
        filename = bgenfiles[ic],
        offsets  = info$file_start_position,
        BM       = G,
        ind_row  = ind_row - 1L,
        ind_col  = ind.col,
        decode   = as.raw(207 - round(0:510 * 100 / 255)),
        dosage   = dosage,
        N        = N,
        ncores   = ncores
      )

      # Return variant info
      info %>%
        dplyr::bind_cols(varinfo) %>%
        dplyr::select(chromosome, marker.ID = ID, rsid, physical.pos = position,
                      allele1, allele2, freq = FREQ, info = INFO)
    }))
  })

  # Create the bigSNP object
  snp.list <- structure(list(genotypes = G, map = snp.info), class = "bigSNP")

  # save it and return the path of the saved object
  rds <- sub_bk(G$backingfile, ".rds")
  saveRDS(snp.list, rds)
  rds
}

################################################################################
