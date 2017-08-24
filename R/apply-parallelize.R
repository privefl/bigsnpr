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
#' @param FUN The function to be applied. It must take a
#' [FBM.code256][FBM.code256-class] as first argument and `ind.chr`,
#' an another argument to provide subsetting over SNPs.
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
#' @examples
#' # parallelize over chromosomes made easy
#' # examples of functions from this package
#' snp_pruning
#' snp_clumping
#' snp_fastImpute
#'
snp_split <- function(infos.chr, FUN, combine, ncores = 1, ...) {

  check_args()
  assert_args(FUN, "ind.chr")

  ind.chrs <- split(seq_along(infos.chr), infos.chr)
  ord.chrs <- order(sapply(ind.chrs, length), decreasing = TRUE)
  inv.ord.chrs <- match(seq_along(ord.chrs), ord.chrs)

  if (ncores == 1) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  # apply the function in decreasing order of chromosome lengths
  res.noorder <- foreach(ic = ord.chrs) %dopar% {
    ind.chr <- ind.chrs[[ic]]
    attr(ind.chr, "chr") <- as.numeric(names(ind.chrs)[ic])

    FUN(ind.chr = ind.chr, ...)
  }
  # reorder the results and combine
  foreach(ic = res.noorder[inv.ord.chrs], .combine = combine) %do% {
    ic
  }
}

################################################################################
