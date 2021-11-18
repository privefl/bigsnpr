################################################################################

ld0 <- function(Gna,
                ind.row = rows_along(Gna),
                ind.col = cols_along(Gna),
                size = 500,
                infos.pos = NULL,
                ncores = 1) {

  if (is.null(infos.pos)) infos.pos <- 1000 * seq_along(ind.col)
  assert_lengths(infos.pos, ind.col)
  assert_sorted(infos.pos)

  ld_scores(
    obj    = Gna,
    rowInd = ind.row,
    colInd = ind.col,
    size   = size * 1000,
    pos    = infos.pos,
    ncores = ncores
  )
}

################################################################################

#' LD scores
#'
#' @inheritParams snp_cor
#'
#' @return A vector of LD scores. For each variant, this is the sum of squared
#'   correlations with the neighboring variants (including itself).
#'
#' @examples
#' test <- snp_attachExtdata()
#' G <- test$genotypes
#'
#' (ld <- snp_ld_scores(G, ind.col = 1:1000))
#'
#' @export
#'
snp_ld_scores <- function(Gna,
                          ind.row = rows_along(Gna),
                          ind.col = cols_along(Gna),
                          size = 500,
                          infos.pos = NULL,
                          ncores = 1) {

  args <- as.list(environment())

  check_args()

  do.call(ld0, args)
}

################################################################################

#' @rdname snp_ld_scores
#' @export
bed_ld_scores <- function(obj.bed,
                          ind.row = rows_along(obj.bed),
                          ind.col = cols_along(obj.bed),
                          size = 500,
                          infos.pos = NULL,
                          ncores = 1) {

  args <- as.list(environment())
  names(args)[names(args) == "obj.bed"] <- "Gna"

  check_args()

  do.call(ld0, args)
}

################################################################################
