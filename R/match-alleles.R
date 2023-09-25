################################################################################

flip_strand <- function(allele) {

  assert_package("dplyr")

  dplyr::case_when(
    allele == "A" ~ "T",
    allele == "C" ~ "G",
    allele == "T" ~ "A",
    allele == "G" ~ "C",
    TRUE ~ NA_character_
  )
}

#' Match alleles
#'
#' Match alleles between summary statistics and SNP information.
#' Match by ("chr", "a0", "a1") and ("pos" or "rsid"), accounting for possible
#' strand flips and reverse reference alleles (opposite effects).
#'
#' @param sumstats A data frame with columns "chr", "pos", "a0", "a1" and "beta".
#' @param info_snp A data frame with columns "chr", "pos", "a0" and "a1".
#' @param strand_flip Whether to try to flip strand? (default is `TRUE`)
#'   If so, ambiguous alleles A/T and C/G are removed.
#' @param join_by_pos Whether to join by chromosome and position (default),
#'   or instead by rsid.
#' @param remove_dups Whether to remove duplicates (same physical position)?
#'   Default is `TRUE`.
#' @param match.min.prop Minimum proportion of variants in the smallest data
#'   to be matched, otherwise stops with an error. Default is `20%`.
#' @param return_flip_and_rev Whether to return internal boolean variables
#'   `"_FLIP_"` (whether the alleles must be flipped: A <--> T & C <--> G,
#'   because on the opposite strand) and `"_REV_"` (whether alleles must be
#'   swapped: `$a0` <--> `$a1`, in which case corresponding `$beta` are multiplied
#'   by -1). Default is `FALSE`.
#'
#' @return A single data frame with matched variants. Values in column `$beta`
#'   are multiplied by -1 for variants with alleles reversed (i.e. swapped).
#'   New variable `"_NUM_ID_.ss"` returns the corresponding row indices of the
#'   input `sumstats` (first argument of this function), and `"_NUM_ID_"`
#'   corresponding to the input `info_snp` (second argument).
#' @export
#'
#' @seealso [snp_modifyBuild]
#'
#' @import data.table
#'
#' @example examples/example-match.R
snp_match <- function(sumstats, info_snp,
                      strand_flip = TRUE,
                      join_by_pos = TRUE,
                      remove_dups = TRUE,
                      match.min.prop = 0.2,
                      return_flip_and_rev = FALSE) {

  sumstats <- as.data.frame(sumstats)
  info_snp <- as.data.frame(info_snp)

  sumstats$`_NUM_ID_` <- rows_along(sumstats)
  info_snp$`_NUM_ID_` <- rows_along(info_snp)

  min_match <- match.min.prop * min(nrow(sumstats), nrow(info_snp))

  join_by <- c("chr", NA, "a0", "a1")
  join_by[2] <- `if`(join_by_pos, "pos", "rsid")

  if (!all(c(join_by, "beta") %in% names(sumstats)))
    stop2("Please use proper names for variables in 'sumstats'. Expected '%s'.",
          paste(c(join_by, "beta"), collapse = ", "))
  if (!all(c(join_by, "pos") %in% names(info_snp)))
    stop2("Please use proper names for variables in 'info_snp'. Expected '%s'.",
          paste(unique(c(join_by, "pos")), collapse = ", "))

  message2("%s variants to be matched.", format(nrow(sumstats), big.mark = ","))

  # first filter to fasten
  sumstats <- sumstats[vctrs::vec_in(sumstats[, join_by[1:2]],
                                     info_snp[, join_by[1:2]]), ]
  if (nrow(sumstats) == 0)
    stop2("No variant has been matched.")

  # augment dataset to match reverse alleles
  if (strand_flip) {
    is_ambiguous <- with(sumstats, paste(a0, a1) %in% c("A T", "T A", "C G", "G C"))
    message2("%s ambiguous SNPs have been removed.",
             format(sum(is_ambiguous), big.mark = ","))
    sumstats2 <- sumstats[!is_ambiguous, ]
    sumstats3 <- sumstats2
    sumstats2$`_FLIP_` <- FALSE
    sumstats3$`_FLIP_` <- TRUE
    sumstats3$a0 <- flip_strand(sumstats2$a0)
    sumstats3$a1 <- flip_strand(sumstats2$a1)
    sumstats3 <- rbind(sumstats2, sumstats3)
  } else {
    sumstats3 <- sumstats
    sumstats3$`_FLIP_` <- FALSE
  }

  sumstats4 <- sumstats3
  sumstats3$`_REV_` <- FALSE
  sumstats4$`_REV_` <- TRUE
  sumstats4$a0 <- sumstats3$a1
  sumstats4$a1 <- sumstats3$a0
  sumstats4$beta <- -sumstats3$beta
  sumstats4 <- rbind(sumstats3, sumstats4)

  matched <- merge(as.data.table(sumstats4), as.data.table(info_snp),
                   by = join_by, all = FALSE, suffixes = c(".ss", ""))

  if (remove_dups) {
    dups <- vctrs::vec_duplicate_detect(matched[, c("chr", "pos")])
    if (any(dups)) {
      matched <- matched[!dups, ]
      message2("Some duplicates were removed.")
    }
  }

  message2("%s variants have been matched; %s were flipped and %s were reversed.",
           format(nrow(matched),         big.mark = ","),
           format(sum(matched$`_FLIP_`), big.mark = ","),
           format(sum(matched$`_REV_`),  big.mark = ","))

  if (nrow(matched) < min_match)
    stop2("Not enough variants have been matched.")

  if (!return_flip_and_rev) matched <- matched[, c("_FLIP_", "_REV_") := NULL]

  as.data.frame(matched[order(chr, pos)])
}

################################################################################

#' Determine reference divergence
#'
#' Determine reference divergence while accounting for strand flips.
#' **This does not remove ambiguous alleles.**
#'
#' @param ref1 The reference alleles of the first dataset.
#' @param alt1 The alternative alleles of the first dataset.
#' @param ref2 The reference alleles of the second dataset.
#' @param alt2 The alternative alleles of the second dataset.
#'
#' @return A logical vector whether the references alleles are the same.
#'   Missing values can result from missing values in the inputs or from
#'   ambiguous matching (e.g. matching A/C and A/G).
#' @export
#'
#' @seealso [snp_match()]
#'
#' @examples
#' same_ref(ref1 = c("A", "C", "T", "G", NA),
#'          alt1 = c("C", "T", "C", "A", "A"),
#'          ref2 = c("A", "C", "A", "A", "C"),
#'          alt2 = c("C", "G", "G", "G", "A"))
same_ref <- function(ref1, alt1, ref2, alt2) {

  # ACTG <- c("A", "C", "T", "G")
  # REV_ACTG <- stats::setNames(c("T", "G", "A", "C"), ACTG)
  #
  # decoder <- expand.grid(list(ACTG, ACTG, ACTG, ACTG)) %>%
  #   dplyr::mutate(status = dplyr::case_when(
  #     # BAD: same reference/alternative alleles in a dataset
  #     (Var1 == Var2) | (Var3 == Var4) ~ NA,
  #     # GOOD/TRUE: same reference/alternative alleles between datasets
  #     (Var1 == Var3) & (Var2 == Var4) ~ TRUE,
  #     # GOOD/FALSE: reverse reference/alternative alleles
  #     (Var1 == Var4) & (Var2 == Var3) ~ FALSE,
  #     # GOOD/TRUE: same reference/alternative alleles after strand flip
  #     (REV_ACTG[Var1] == Var3) & (REV_ACTG[Var2] == Var4) ~ TRUE,
  #     # GOOD/FALSE: reverse reference/alternative alleles after strand flip
  #     (REV_ACTG[Var1] == Var4) & (REV_ACTG[Var2] == Var3) ~ FALSE,
  #     # BAD: the rest
  #     TRUE ~ NA
  #   )) %>%
  #   reshape2::acast(Var1 ~ Var2 ~ Var3 ~ Var4, value.var = "status")

  decoder <- structure(
    c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE,
      NA, NA, FALSE, NA, NA, NA, NA, NA, NA, TRUE, NA, NA, FALSE, NA, NA, NA,
      TRUE, NA, NA, NA, NA, NA, FALSE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
      TRUE, NA, NA, FALSE, NA, NA, TRUE, NA, NA, FALSE, NA, NA, NA, NA, FALSE,
      NA, NA, TRUE, NA, NA, NA, NA, NA, NA, FALSE, NA, NA, TRUE, NA, NA, NA, NA,
      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, FALSE, NA,
      NA, TRUE, NA, NA, FALSE, NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA,
      NA, TRUE, NA, NA, NA, NA, NA, FALSE, NA, NA, NA, NA, FALSE, NA, NA, NA,
      NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, NA, NA, FALSE,
      NA, NA, TRUE, NA, NA, FALSE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
      NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, NA, NA, FALSE, NA, NA, NA, NA,
      NA, NA, TRUE, NA, NA, FALSE, NA, NA, NA, NA, FALSE, NA, NA, TRUE, NA, NA,
      FALSE, NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, FALSE, NA,
      NA, NA, NA, NA, TRUE, NA, NA, NA, FALSE, NA, NA, TRUE, NA, NA, NA, NA, NA,
      NA, FALSE, NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
      NA, NA, NA, NA, NA),
    .Dim = rep(4, 4), .Dimnames = rep(list(c("A", "C", "T", "G")), 4)
  )

  to_decode <- do.call("cbind", lapply(list(ref1, alt1, ref2, alt2), as.character))
  decoder[to_decode]
}

################################################################################
