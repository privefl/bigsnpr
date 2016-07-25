#' @title Modify a "bigSNP" and save these modifications.
#' @description Modify a \code{bigSNP} and save these modifications
#' in the corresponding .rds file.
#' @param x A \code{bigSNP}.
#' @param to.save Is the modification to be saved ?
#' @return Nothing.
#' @name modif-save



#' @rdname modif-save
#' @export
SaveModifs <- function(x) {
  saveRDS(x, file.path(x$backingpath, paste0(x$backingfile, ".rds")))

  return()
}


#' @rdname modif-save
#' @param pop.files Character vector of file names where to
#' find the population of the individuals of the study.
#' @param col.sample.ID Index of the column containing the
#' sample IDs to match with those of the study.
#' @param col.family.ID Index of the column containt the populations.
#' @export
GetPops <- function(x,
                    pop.files,
                    col.sample.ID,
                    col.family.ID,
                    to.save = FALSE) {
  obj <- foreach::foreach(f = pop.files, .combine = 'rbind')
  expr_fun <- function(f) data.table::fread(f, data.table = FALSE)
  data.pop <- foreach2(obj, expr_fun, 1)

  num <- match(x$fam$sample.ID, data.pop[, col.sample.ID])
  no.match <- sum(is.na(num))
  if (no.match == 0) {
    printf("Each individual has been matched\n")
  } else {
    printf("There are %d individuals which have not been matched\n")
  }

  x$fam$family.ID = data.pop[num, col.family.ID]

  if (to.save) SaveModifs(x)

  return()
}
