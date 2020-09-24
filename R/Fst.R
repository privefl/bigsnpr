################################################################################

sum_transfo_af <- function(list_df_af, expr) {
  env <- parent.frame()
  Reduce(`+`,
         lapply(list_df_af,
                as.function(c(alist(df_af = ), substitute(expr)), envir = env)))
}

################################################################################

#' Fixation index (Fst)
#'
#' Fixation index (Fst), either per variant, or genome-wide
#'
#' @param list_df_af List of data frames with `$af` (allele frequency per variant)
#'   and `$N` (sample size per variant). Typically, the outputs of [bed_MAF()].
#'   Each new data frame of the list should correspond to a different population.
#' @param min_maf Minimum MAF threshold (for the average of populations) to be
#'   included in the final results. Default is `0` (remove monomorphic variants).
#' @param overall Whether to compute Fst genome-wide (`TRUE`) or per variant
#'   (`FALSE`, the default).
#'
#' @return If `overall`, then one value, otherwise a value for each variant
#'   with missing values for the variants not passing `min_maf`.
#'   This should be equivalent to using '`--fst --within`' in PLINK.
#'
#' @export
#'
#' @references
#' Weir, B. S., & Cockerham, C. C. (1984). Estimating F-statistics for the
#' analysis of population structure. Evolution, 1358-1370.
#'
#' @examples
#' bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#' obj.bed <- bed(bedfile)
#'
#' pop <- rep(1:3, c(143, 167, 207))
#' ind_pop <- split(seq_along(pop), pop)
#' list_df_af <- lapply(ind_pop, function(ind) bed_MAF(obj.bed, ind.row = ind))
#'
#' snp_fst(list_df_af)
#' snp_fst(list_df_af[c(1, 2)], overall = TRUE)
#' snp_fst(list_df_af[c(1, 3)], overall = TRUE)
#' snp_fst(list_df_af[c(3, 2)], overall = TRUE)
#'
snp_fst <- function(list_df_af, min_maf = 0, overall = FALSE) {

  r <- length(list_df_af)
  if (r < 2)
    stop2("You should provide frequencies for at least 2 populations.")
  if (min_maf < 0 | min_maf > 0.45)
    stop2("Parameter 'min_maf' should be in range [0, 0.45].")
  lapply(list_df_af, function(df_af)
    assert_df_with_names(df_af, c("af", "N")))

  n_sum <- sum_transfo_af(list_df_af, df_af$N)
  n_bar <- n_sum / r
  n_sqsum <- sum_transfo_af(list_df_af, df_af$N^2)
  n_c <- (n_sum - n_sqsum / n_sum) / (r - 1)

  af_n_sum <- sum_transfo_af(list_df_af, df_af$af * df_af$N)
  p_bar <- af_n_sum / n_sum

  diff_af_n_sum <- sum_transfo_af(list_df_af, (df_af$af - p_bar)^2 * df_af$N)
  s2 <- diff_af_n_sum / n_bar / (r - 1)

  h_n_sum <- sum_transfo_af(list_df_af, 2 * df_af$af * (1 - df_af$af) * df_af$N)
  h_bar <- h_n_sum / n_sum

  a <- n_bar / n_c * (s2 - 1 / (n_bar - 1) *
                        (p_bar * (1 - p_bar) - (r - 1) / r * s2 - h_bar / 4))
  b <- n_bar / (n_bar - 1) *
    (p_bar * (1 - p_bar) - (r - 1) / r * s2 - (2 * n_bar - 1) / (4 * n_bar) * h_bar)
  c <- h_bar / 2

  keep <- (p_bar > min_maf & p_bar < (1 - min_maf))

  if (overall) {
    ind_keep <- which(keep)
    sum(a[ind_keep]) / sum((a + b + c)[ind_keep])
  } else {
    ifelse(keep, a / (a + b + c), NA)
  }
}

################################################################################
