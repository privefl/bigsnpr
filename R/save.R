################################################################################

#' @title Modify a "bigSNP" and save these modifications.
#' @description Modify a \code{bigSNP} and save these modifications
#' in the corresponding .rds file. As \code{bigSNP} is an S3 class, you
#' can add any slot you want to this class, then use \code{SaveModifs} to
#' save modifications in the corresponding .rds backing file.
#' @param x A \code{bigSNP}.
#' @return The (modified) \code{bigSNP}.
#' @name modif-save
#' @example examples/example.save.R
#' @seealso \code{\link{rle}}
NULL

################################################################################


#' @rdname modif-save
#' @export
SaveModifs <- function(x) {
  saveRDS(x, file.path(x$backingpath, paste0(x$backingfile, ".rds")))

  return(x)
}

################################################################################

#' @rdname modif-save
#' @description \code{GetPops}: Get populations of individuals by matching
#' from an external file or get them directly from the slot "family.ID" of
#' the slot "fam".
#' It will add an slot called "pop" to the slot "fam".
#' @param pop.files Character vector of file names where to
#' find the population of the individuals of the study.
#' Default just copy the slot "family.ID" to the slot "pop".
#' @param col.sample.ID Index of the column containing the
#' sample IDs to match with those of the study.
#' @param col.family.ID Index of the column containt the populations.
#' @param ... Any additional parameter to pass to \code{fread}.
#' @seealso \code{\link[data.table]{fread}}
#' @export
GetPops <- function(x,
                    pop.files = NULL,
                    col.sample.ID,
                    col.family.ID,
                    ...) {
  if (is.null(pop.files)) {
    x$fam$pop <- x$fam$family.ID
  } else {
    obj <- foreach::foreach(f = pop.files, .combine = 'rbind')
    expr_fun <- function(f) data.table::fread(f, data.table = FALSE)
    data.pop <- foreach2(obj, expr_fun, 1)

    num <- match(x$fam$sample.ID, data.pop[, col.sample.ID])
    no.match <- sum(is.na(num))
    if (no.match != 0) {
      warning(sprintf(
        "There are %d individuals which have not been matched", no.match))
    }

    x$fam$pop <- data.pop[num, col.family.ID]
  }

  SaveModifs(x)
}

################################################################################

#' @rdname modif-save
#' @description \code{GetPhenos}: Modify phenotypes to be -1 (unaffected),
#' 1 (affected) or NA (missing).
#' It will add an slot called "pheno" to the slot "fam".
#' @param coded01 Are the phenotypes coded 0 (unaffected) / 1 (affected)
#' rather than respectively 1 and 2? Default is \code{FALSE}.
#' @export
GetPhenos <- function(x, coded01 = FALSE) {
  tmp <- x$fam$affection

  if (coded01) {
    noNAs <- tmp %in% c(0, 1)
    tmp[!noNAs] <- NA
    x$fam$pheno <- tmp * 2 - 1
  } else {
    noNAs <- tmp %in% c(1, 2)
    tmp[!noNAs] <- NA
    x$fam$pheno <- tmp * 2 - 3
  }

  no.match <- sum(!noNAs)
  if (no.match != 0)
    warning(sprintf("There are %d missing phenotypes", no.match))

  SaveModifs(x)
}

################################################################################
