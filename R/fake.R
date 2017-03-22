################################################################################

#' Fake a "bigSNP"
#'
#' @param n Number of individuals.
#' @param m Number of SNPs.
#'
#' @return A new temporary `bigSNP` object representing `n` individuals and
#' `m` SNPs. The genotype `big.matrix` is initialized with missing values.
#'
#' @keywords internal
#'
#' @examples
#' N <- 5
#' M <- 12
#' test <- snp_fake(N, M)
#' str(test)
#'
#' # The genotype `big.matrix` is initialized with missing values
#' X <- attach.BM(test$genotypes)
#' X[,]
#'
#' # Modify the genotype `big.matrix`
#' X[] <- sample(as.raw(0:3), size = length(X), replace = TRUE)
#' X[,]
#'
#' rm(X)
#'
#' @export
snp_fake <- function(n, m) {

  # backingfiles
  tmpfile <- tempfile()
  backingfile <- basename(tmpfile)
  backingpath <- dirname(tmpfile)

  # constructing a fake genotype big.matrix
  bigGeno <- big.matrix(n, m, type = "raw", init = as.raw(3),
                  backingfile = paste0(backingfile, ".bk"),
                  backingpath = backingpath,
                  descriptorfile = paste0(backingfile, ".desc"))
  bigGeno.code <- as.BM.code(bigGeno, code = CODE_012)

  # fam
  fam <- data.frame(0L, paste0("ind_", 1:n), 0L, 0L, 0L, -9L,
                    stringsAsFactors = FALSE)
  names(fam) <- NAMES.FAM

  # map
  map <- data.frame(1L, paste0("snp_", 1:m), 0L, 0L,
                    ifelse(cond <- (runif(m) > 0.5), "A", "T"),
                    ifelse(!cond, "A", "T"),
                    stringsAsFactors = FALSE)
  names(map) <- NAMES.MAP

  # create the `bigSNP`, save it and return it
  rds <- paste0(tmpfile, ".rds")
  snp_list <- structure(list(genotypes = describe(bigGeno.code),
                             fam = fam,
                             map = map,
                             savedIn = rds),
                        class = "bigSNP")
  saveRDS(snp_list, rds)
  snp_list
}

################################################################################
