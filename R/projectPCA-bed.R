################################################################################

#' Download 1000G
#'
#' Download 1000 genomes project (phase 3) data in PLINK bed/bim/fam format,
#' including 2493 (mostly unrelated) individuals
#' and ~1.4M SNPs in common with HapMap3.
#'
#' @param dir The directory where to put the downloaded files.
#' @param overwrite Whether to overwrite files when downloading and unzipping?
#'   Default is `FALSE`, so this function does not do anything the second time.
#'
#' @return The path of the downloaded bed file.
#'
#' @export
#'
download_1000G <- function(dir, overwrite = FALSE) {

  zip <- file.path(dir, "1000G_phase3_common_hapmap.zip")
  if (overwrite || !file.exists(zip)) {
    utils::download.file("https://ndownloader.figshare.com/files/16775654",
                         destfile = zip)
  }

  files_in_zip <- paste0("1000G_phase3_common_hapmap_norel",
                         c(".bed", ".bim", ".fam", ".fam2"))
  files_unzipped <- file.path(dir, files_in_zip)

  if (!overwrite) files_in_zip <- files_in_zip[!file.exists(files_unzipped)]

  if (length(files_in_zip) > 0)
    utils::unzip(zip, exdir = dir, files = files_in_zip)

  files_unzipped[1]
}

################################################################################
