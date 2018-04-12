################################################################################

#' Get genes
#'
#' Get genes associated with SNPs.
#'
#' @param rsid A character vector of 'rs' ID of SNPs to investigate.
#' @inheritParams bigsnpr-package
#'
#' @return A character vector of genes in the form `"<name>:<ID>". Note that
#'   there can be multiple genes per SNP (separated by a comma), or none (`NA`).
#' @export
#'
#' @examples
#' rsid <- c("rs3934834", "rs3737728", "rs6687776", "rs9651273", "rs4970405",
#'           "rs12726255", "rs2298217", "rs4970362", "rs9660710", "rs4970420")
#' snp_gene(rsid)
snp_gene <- function(rsid, ncores = 1) {

  if (requireNamespace("rsnps", quietly = TRUE)) {
    big_apply(FBM(0, 0), a.FUN = function(X, ind, ID) {
      rsnps::ncbi_snp_summary(ID[ind])$gene2
    }, a.combine = 'c', ind = seq_along(rsid),
    block.size = 300, ncores = ncores, ID = rsid)
  } else {
    stop2(paste("Please install R package {rsnps}",
                "with `devtools::install_github('ropensci/rsnps')`."))
  }
}

################################################################################
