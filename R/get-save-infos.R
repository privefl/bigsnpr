################################################################################

#' Save modifications
#'
#' Save a \code{bigSNP} after having made some modifications to it.
#' As \code{bigSNP} is an S3 class, you can add any slot you want
#' to an object of this class, then use \code{snp_save} to
#' save these modifications in the corresponding ".rds" backing file.
#'
#' @inheritParams bigsnpr-package
#'
#' @return The (saved) \code{bigSNP}.
#'
#' @example examples/example-save.R
#'
#' @export
snp_save <- function(x) {
  saveRDS(x, sub("\\.bk$", ".rds", x$genotypes$backingfile))

  x
}

################################################################################

#' Get sample information
#'
#' Get information of individuals by matching from an external file.
#'
#' @inheritParams bigsnpr-package
#' @param df.or.files Either
#' - A `data.frame`,
#' - A character vector of file names where to find at the information you want.
#'   You should have one column for family IDs and one for sample IDs.
#' @param col.sample.ID Index of the column containing the
#' sample IDs to match with those of the study. Default is the first one.
#' @param col.family.ID Index of the column containing the
#' family IDs to match with those of the study. Default is the second one.
#' @param col.infos Indices of the column containing the information you want.
#' Default is all but the first and the second columns.
#' @param pair.sep Separator used for concatenation family and sample IDs
#' in order to match easier. Default is `"-_-"`.
#' @param ... Any additional parameter to pass to [fread][data.table::fread].
#' Particularly, option `header = FALSE` is sometimes needed.
#'
#' @return The requested information as a `data.frame`.
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
#' infos <- snp_getSampleInfos(test, files, header = FALSE)
#' rle(infos[[1]])
#'
#' @seealso [list.files]
#'
#' @export
snp_getSampleInfos <- function(x, df.or.files,
                               col.family.ID = 1,
                               col.sample.ID = 2,
                               col.infos = -c(1, 2),
                               pair.sep = "-_-",
                               ...) {

  check_args()

  if (is.data.frame(df.or.files)) {
    data.infos <- df.or.files
  } else if (is.character(df.or.files)) {
    data.infos <- foreach(f = df.or.files, .combine = 'rbind') %do% {
      data.table::fread(f, data.table = FALSE, ...)
    }
  } else {
    stop2("'df.or.files' must be a data.frame or a vector of file paths.")
  }

  to.match   <- paste(x$fam$family.ID, x$fam$sample.ID, sep = pair.sep)
  from.match <- paste(data.infos[, col.family.ID],
                      data.infos[, col.sample.ID],
                      sep = pair.sep)
  num <- match(to.match, from.match)
  if (no.match <- sum(is.na(num)))
    warning2("There are %d individuals which have not been matched", no.match)

  data.infos[num, col.infos, drop = FALSE]
}

################################################################################
