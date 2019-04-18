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
#'   ambiguous strands.
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

  decoder[cbind(ref1, alt1, ref2, alt2)]
}

################################################################################

flip_strand <- function(allele) {

  if (!requireNamespace("dplyr", quietly = TRUE))
    stop2("Please install package 'dplyr'.")

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
#' Match by "chr", "pos", "a0", "a1", accounting for strand flips (as an option)
#' and reverse reference alleles (opposite effect).
#'
#' @param sumstats A data frame with columns "chr", "pos", "a0", "a1" and "beta".
#' @param info_snp A data frame with columns "chr", "pos", "a0" and "a1".
#' @param strand_flip Whether it should try to flip strand? (default) If so,
#'  ambiguous alleles A/T and C/G are removed.
#'
#' @return A single data frame joined by chromosome, position and alleles.
#' @export
#'
#' @import data.table
#'
#' @example examples/example-match.R
snp_match <- function(sumstats, info_snp, strand_flip = TRUE) {

  if (!all(c("chr", "pos", "a0", "a1", "beta") %in% names(sumstats)))
    stop2("Please use proper names for variables in 'sumstats'.")
  if (!all(c("chr", "pos", "a0", "a1") %in% names(info_snp)))
    stop2("Please use proper names for variables in 'info_snp'.")

  message2("%s variants in summary statistics.",
           format(nrow(sumstats), big.mark = ","))

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
                   by = c("chr", "pos", "a0", "a1"), all = FALSE)
  message2("%s variants have been matched; %s were flipped and %s were reversed.",
           format(nrow(matched),         big.mark = ","),
           format(sum(matched$`_FLIP_`), big.mark = ","),
           format(sum(matched$`_REV_`),  big.mark = ","))

  as.data.frame(matched[, c("_FLIP_", "_REV_") := NULL][order(chr, pos)])
}

################################################################################
