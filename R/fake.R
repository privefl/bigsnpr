################################################################################

#' Fake a "bigSNP"
#'
#' @param n Number of individuals.
#' @param m Number of SNPs.
#'
#' @return A new temporary `bigSNP` object representing `n` individuals
#'   and `m` SNPs. The genotype Filebacked Big Matrix is initialized
#'   with missing values.
#'
#' @keywords internal
#'
#' @examples
#' (test <- snp_fake(5, 12))
#'
#' # The genotype Filebackeg Big Matrix is initialized with missing values
#' G <- test$genotypes
#' G[]
#'
#' # Modify the genotype `big.matrix`
#' G[] <- sample(as.raw(0:3), size = length(G), replace = TRUE)
#' G[]
#'
#' @export
#'
snp_fake <- function(n, m) {

  # constructing a fake genotype big.matrix
  bigGeno <- FBM.code256(n, m, code = CODE_012, init = as.raw(3))

  # fam
  fam <- data.frame(0L, paste0("ind_", 1:n), 0L, 0L, 0L, -9L,
                    stringsAsFactors = FALSE)
  names(fam) <- NAMES.FAM

  # map
  map <- data.frame(1L, paste0("snp_", 1:m), 0L, 0L,
                    ifelse(cond <- (stats::runif(m) > 0.5), "A", "T"),
                    ifelse(!cond, "A", "T"),
                    stringsAsFactors = FALSE)
  names(map) <- NAMES.MAP

  # Return the bigSNP object
  structure(list(genotypes = bigGeno,
                 fam = fam,
                 map = map),
            class = "bigSNP")
}

################################################################################
