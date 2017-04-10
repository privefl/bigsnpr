################################################################################

#' Save modifications
#'
#' Save a \code{bigSNP} after having made some modifications to it.
#' As \code{bigSNP} is an S3 class, you can add any slot you want
#' to an object of this class, then use \code{snp_saveModifs} to
#' save these modifications in the corresponding ".rds" backing file.
#'
#' @inheritParams bigsnpr-package
#'
#' @return The (saved) \code{bigSNP}.
#'
#' @example examples/example-saveModifs.R
#'
#' @export
snp_saveModifs <- function(x) {
  saveRDS(x, x$savedIn)

  x
}

################################################################################

#' Get sample infos
#'
#' Get information of individuals by matching from an external file.
#'
#' @param files Character vector of file names where to find at the information
#' you want. You should have one column for family IDs and one for sample IDs.
#' @param col.sample.ID Index of the column containing the
#' sample IDs to match with those of the study.
#' @param col.family.ID Index of the column containing the
#' family IDs to match with those of the study.
#' @param col.infos Indices of the column containing the information you want.
#' @param pair.sep Separator used for concataining family and sample IDs
#' in order to match easier. Default is `"-_-"`.
#' @param ... Any additional parameter to pass to [fread][data.table::fread].
#' Particularly, option `header = FALSE` is sometimes needed.
#'
#' @return The requested information. It will be a vector if you request only
#' one column, and a `data.frame` otherwise.
#' @import foreach
#'
#' @examples
#' test <- snp_attachExtdata()
#' # Just after reading
#' rle(test$fam$family.ID)
#' # Get populations clusters from external files
#' files <- system.file("extdata", paste0("cluster", 1:3), package = "bigsnpr")
#' data.table::fread(files[1])
#' # need header option
#' data.table::fread(files[1], header = FALSE)
#' infos <- snp_getSampleInfos(test, files = files,
#'                             col.family.ID = 1,
#'                             col.sample.ID = 2,
#'                             col.infos = 3,
#'                             header = FALSE)
#' rle(infos)
#'
#' @export
snp_getSampleInfos <- function(x, files,
                               col.family.ID,
                               col.sample.ID,
                               col.infos,
                               pair.sep = "-_-",
                               ...) {

  data.infos <- foreach(f = files, .combine = 'rbind') %do% {
    data.table::fread(f, data.table = FALSE, ...)
  }

  to.match   <- paste(x$fam$family.ID, x$fam$sample.ID, sep = pair.sep)
  from.match <- paste(data.infos[, col.family.ID],
                      data.infos[, col.sample.ID],
                      sep = pair.sep)
  num <- match(to.match, from.match)
  if (no.match <- sum(is.na(num)))
    warning2("There are %d individuals which have not been matched", no.match)

  data.infos[num, col.infos]
}

################################################################################
