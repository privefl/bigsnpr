################################################################################

#' Title
#'
#' @param n
#' @param m
#'
#' @return
#' @export
#'
#' @examples
snp_fake <- function(n, m) {
  # backingfiles
  tmpfile <- tempfile()
  backingfile <- basename(tmpfile)
  backingpath <- dirname(tmpfile)

  # constructing a fake genotype big.matrix
  X <- big.matrix(n, m, type = "char", init = NA,
                  backingfile = paste0(backingfile, ".bk"),
                  backingpath = backingpath,
                  descriptorfile = paste0(backingfile, ".desc"))

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

  structure(list(genotypes = X,
                 fam = fam,
                 map = map,
                 backingfile = backingfile,
                 backingpath = backingpath),
            class = "bigSNP")
}

################################################################################
