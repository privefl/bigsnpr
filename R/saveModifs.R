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

#' Get population infos
#'
#' Get populations of individuals by matching from an external file and
#' __overwrite__ the slot __`family.ID`__ of the slot `fam`.
#'
#' @param pop.files Character vector of file names where to
#' find the population of the individuals of the study.
#' @param col.sample.ID Index of the column containing the
#' sample IDs to match with those of the study.
#' @param col.family.ID Index of the column containt the populations.
#' @param save Save modifications in the ".rds" file? Default is `FALSE`.
#' @param ... Any additional parameter to pass to [fread][data.table::fread].
#'
#' @return The modified `bigSNP`. If `save`, modifications are saved
#' in the ".rds" file.
#' @import foreach
#'
#' @examples
#' test <- snp_attachExtdata()
#'
#' # Just after reading
#' print(rle(test$fam$family.ID))
#'
#' # Get populations clusters from external files
#' files <- system.file("extdata", paste0("cluster", 1:3), package = "bigsnpr")
#' data.table::fread(files[1])
#' test <- snp_getPops(test, pop.files = files,
#'                     col.sample.ID = 2,
#'                     col.family.ID = 3)
#'
#' rle(test$fam$family.ID)
#'
#' @export
snp_getPops <- function(x,
                        pop.files,
                        col.sample.ID,
                        col.family.ID,
                        save = FALSE,
                        ...) {
  data.pop <- foreach(f = pop.files, .combine = 'rbind') %do% {
    data.table::fread(f, data.table = FALSE, ...)
  }

  num <- match(x$fam$sample.ID, data.pop[, col.sample.ID])
  if (no.match <- sum(is.na(num)))
    warning2("There are %d individuals which have not been matched", no.match)

  x$fam$family.ID <- data.pop[num, col.family.ID]

  `if`(save, snp_saveModifs(x), x)
}

################################################################################
