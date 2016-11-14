source("R/utils.R")

requestSNPs <- function(snps) {
  m <- length(snps)
  intervals <- CutBySize(m, 500)
  obj <- foreach::foreach(i = 1:nrow(intervals), 
                          .combine = 'rbind')
  expr_fun <- function(i) {
    rsnps::NCBI_snp_query(snps[seq2(intervals[i, ])])
  }
  foreach2(obj, expr_fun, ncores = 1)
}

infos <- requestSNPs(snps)