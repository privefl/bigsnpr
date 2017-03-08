################################################################################

#' Split-parApply-Combine
#'
#' A Split-Apply-Combine strategy to parallelize the evaluation of a function
#' on each SNP, independently.
#'
#' This function splits indices for each chromosome, then apply a given function
#' to each part (chromosome) and finally combine the results.
#'
#' @inheritParams bigsnpr-package
#' @param infos.chr Vector of integers specifying the chromosome of each SNP.
#' Typically the slot `chromosome` of the slot `map` of a `bigSNP`.
#' @param FUN The function to be applied. It must take a
#' [BM.code][BM.code-class] as first argument and `ind.chr`, an another argument
#' to provide subsetting over SNPs.
#' You can access the number of the chromosome by using `attr(ind.chr, "chr")`.
#' @param combine function that is used by [foreach] to process the tasks
#' results as they generated. This can be specified as either a function or a
#' non-empty character string naming the function. Specifying 'c' is useful
#' for concatenating the results into a vector, for example. The values 'cbind'
#' and 'rbind' can combine vectors into a matrix. The values '+' and '*' can be
#' used to process numeric data. By default, the results are returned in a list.
#' @param ... Extra arguments to be passed to `FUN`.
#'
#' @return The result of [foreach].
#' @export
#' @import foreach
#'
#' @example
snp_split <- function(X.desc, infos.chr, FUN, combine, ncores = 1, ...) {

  ind.chrs <- split(seq_along(infos.chr), infos.chr)

  if (ncores == 1) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  foreach(ic = seq_along(ind.chrs), .combine = combine) %dopar% {
    ind.chr <- ind.chrs[[ic]]
    attr(ind.chr, "chr") <- as.numeric(names(ind.chrs)[ic])

    FUN(X.desc, ind.chr = ind.chr, ...)
  }
}

################################################################################
