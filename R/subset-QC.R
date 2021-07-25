################################################################################

#' Subset a bigSNP
#'
#' Subset (copy) of a `bigSNP`, also stored on disk.
#'
#' @inheritParams bigsnpr-package
#' @param ind.row Indices of the rows (individuals) to keep.
#'   Negative indices __can__ be used to exclude row indices.
#'   Default: keep them all.
#' @param ind.col Indices of the columns (SNPs) to keep.
#'   Negative indices __can__ be used to exclude column indices.
#'   Default: keep them all.
#' @param backingfile Prefix of the two new files created (".bk" and ".rds").
#'   By default, it is automatically determined by appending "_sub" and a number
#'   to the prefix of the input bigSNP backing files.
#'
#' @export
#' @return The path to the RDS file that stores the `bigSNP` object.
#'
#' @seealso [bigSNP][bigSNP-class]
#' @examples
#' str(test <- snp_attachExtdata())
#'
#' # keep only first 50 samples and SNPs
#' rdsfile <- snp_subset(test, ind.row = 1:50, ind.col = 1:50)
#' str(snp_attach(rdsfile))
#'
#' # remove only first 50 samples and SNPs
#' rdsfile2 <- snp_subset(test, ind.row = -(1:50), ind.col = -(1:50))
#' str(snp_attach(rdsfile2))
#'
snp_subset <- function(x,
                       ind.row = rows_along(x$fam),
                       ind.col = rows_along(x$map),
                       backingfile = NULL) {


  # Support for negative indices
  ind.row <- rows_along(x$fam)[ind.row]
  ind.col <- rows_along(x$map)[ind.col]

  check_args()

  # https://stackoverflow.com/q/19565621/6103040
  new_fam <- x$fam[ind.row, , drop = FALSE]
  rownames(new_fam) <- rows_along(new_fam)
  new_map <- x$map[ind.col, , drop = FALSE]
  rownames(new_map) <- rows_along(new_map)

  if (is.null(backingfile)) backingfile <- getNewFile(x, "sub")

  G <- x$genotypes
  # Create new FBM and fill it
  G2 <- FBM.code256(
    nrow = length(ind.row),
    ncol = length(ind.col),
    code = G$code256,
    init = NULL,
    backingfile = backingfile,
    create_bk = TRUE
  )
  replaceSNP(G2, G, rowInd = ind.row, colInd = ind.col)

  # Create the bigSNP object
  snp.list <- structure(list(genotypes = G2,
                             fam = new_fam,
                             map = new_map),
                        class = "bigSNP")

  # save it and return the path of the saved object
  rds <- sub_bk(G2$backingfile, ".rds")
  saveRDS(snp.list, rds)
  rds
}

################################################################################

#' @export
#' @rdname snp_subset
#' @param ... Not used.
subset.bigSNP <- function(x,
                          ind.row = rows_along(x$fam),
                          ind.col = rows_along(x$map),
                          backingfile = NULL,
                          ...) {

  bigassertr::assert_nodots()

  snp_subset(x, ind.row = ind.row, ind.col = ind.col, backingfile = backingfile)
}

################################################################################
