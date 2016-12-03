################################################################################

#' Modify a "bigSNP"
#'
#' Modify a \code{bigSNP} and save these modifications
#' in the corresponding .rds file. As \code{bigSNP} is an S3 class, you
#' can add any slot you want to this class, then use \code{snp_saveModifs} to
#' save modifications in the corresponding .rds backing file.
#'
#' @inheritParams bigsnpr-package
#' @return The (modified) \code{bigSNP}.
#' @name modif-save
#' @example examples/example.save.R
NULL

################################################################################


#' @rdname modif-save
#' @export
snp_saveModifs <- function(x) {
  saveRDS(x, file.path(x$backingpath, paste0(x$backingfile, ".rds")))

  return(x)
}

################################################################################

#' @rdname modif-save
#' @description `GetPops`: Get populations of individuals by matching
#' from an external file and __overwrite__ the slot __`family.ID`__
#' of the slot `fam`.
#' @param pop.files Character vector of file names where to
#' find the population of the individuals of the study.
#' @param col.sample.ID Index of the column containing the
#' sample IDs to match with those of the study.
#' @param col.family.ID Index of the column containt the populations.
#' @param ... Any additional parameter to pass to [fread][data.table::fread].
#' @export
#' @include foreach
GetPops <- function(x,
                    pop.files = NULL,
                    col.sample.ID,
                    col.family.ID,
                    ...) {
  data.pop <- foreach(f = pop.files, .combine = 'rbind') %do% {
    data.table::fread(f, data.table = FALSE)
  }

  num <- match(x$fam$sample.ID, data.pop[, col.sample.ID])
  no.match <- sum(is.na(num))
  if (no.match != 0) {
    warning(sprintf(
      "There are %d individuals which have not been matched", no.match))
  }

  x$fam$family.ID <- data.pop[num, col.family.ID]

  snp_saveModifs(x)
}

################################################################################

#' #' @rdname modif-save
#' #' @description \code{GetPhenos}: Modify phenotypes to be -1 (unaffected),
#' #' 1 (affected) or NA (missing).
#' #' It will add an slot called "pheno" to the slot "fam".
#' #' @param coded01 Are the phenotypes coded 0 (unaffected) / 1 (affected)
#' #' rather than respectively 1 and 2? Default is \code{FALSE}.
#' #' @export
#' GetPhenos <- function(x, coded01 = FALSE) {
#'   tmp <- x$fam$affection
#'
#'   if (coded01) {
#'     noNAs <- tmp %in% c(0, 1)
#'     tmp[!noNAs] <- NA
#'     x$fam$pheno <- tmp * 2 - 1
#'   } else {
#'     noNAs <- tmp %in% c(1, 2)
#'     tmp[!noNAs] <- NA
#'     x$fam$pheno <- tmp * 2 - 3
#'   }
#'
#'   no.match <- sum(!noNAs)
#'   if (no.match != 0)
#'     warning(sprintf("There are %d missing phenotypes", no.match))
#'
#'   snp_saveModifs(x)
#' }

################################################################################
