
#' Convert matrix to genotypes
#'
#' Convert a standard R matrix to a FBM.code256.
#'
#' @param mat A standard R matrix with
#'
#' @return A temporary FBM.code256.
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' mat <- matrix(NA_real_, 8, 12)
#' mat[] <- sample(c(0, 1, 2), size = length(mat), replace = TRUE)
#' mat
#' G <- mat_to_geno(mat)
#' G
#' G[]
#'
mat_to_geno <- function(mat) {

  suppressWarnings(
    mat %>%
      big_copy(type = "raw") %>%
      add_code256(code = CODE_012)
  )
}

################################################################################
