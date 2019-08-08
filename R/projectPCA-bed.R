################################################################################

#' Download 1000G
#'
#' Download 1000 genomes project (phase 3) data in PLINK bed/bim/fam format,
#' including 2493 (mostly unrelated) individuals
#' and ~1.4M SNPs in common with HapMap3.
#'
#' @param dir The directory where to put the downloaded files.
#' @param overwrite Whether to overwrite files when downloading and unzipping?
#'   Default is `FALSE`.
#' @param delete_zip Whether to delete zip after uncompressing the file in it?
#'   Default is `TRUE`.
#'
#' @return The path of the downloaded bed file.
#'
#' @export
#'
download_1000G <- function(dir, overwrite = FALSE, delete_zip = TRUE) {

  files_unzipped <- file.path(dir, paste0("1000G_phase3_common_norel",
                                          c(".bed", ".bim", ".fam", ".fam2")))

  if (overwrite || !all(file.exists(files_unzipped))) {

    zip <- file.path(dir, "1000G_phase3_common.zip")
    if (overwrite || !file.exists(zip)) {
      utils::download.file("https://ndownloader.figshare.com/files/17016137",
                           destfile = zip)
      if (delete_zip) on.exit(unlink(zip), add = TRUE)
    }
    utils::unzip(zip, exdir = dir)

  }

  files_unzipped[1]
}

################################################################################
