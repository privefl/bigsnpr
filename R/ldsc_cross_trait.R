

################################################################################

#' Cross-Trait LD score regression
#'
#' @param ld_score Vector of LD scores.
#' @param ld_size Number of variants used to compute `ld_score`.
#' @param z1 Vector of z-scores for trait 1.
#' @param z2 Vector of z-scores for trait 2.
#' @param sample_size_1 Sample size of GWAS for trait 1.
#'   Possibly a vector, or just a single value.
#' @param sample_size_2 Sample size of GWAS for trait 2.
#'   Possibly a vector, or just a single value.
#' @param step1_chisq_max Threshold on `chi2` in step 1. Default is `30`.
#' @param chi2_thr2 Threshold on `chi2` in step 2. Default is `Inf` (none).
#' @param blocks Either a single number specifying the number of blocks,
#'   or a vector of integers specifying the block number of each `chi2` value.
#'   Default is `200` for `snp_ldsc()`, dividing into 200 blocks of approximately
#'   equal size. `NULL` can also be used to skip estimating standard errors,
#'   which is the default for `snp_ldsc2()`.
#' @param intercept You can constrain the intercept to some value (e.g. 0).
#'   Default is `NULL` (the intercept is estimated).
#'   Use a value of 0 if you are sure there is no overlap between GWAS samples.
#' @inheritParams bigsnpr-package
#'
#' @return Vector of 4 values (or only the first 2 if `blocks = NULL`):
#'  - `[["int"]]`: LDSC regression intercept,
#'  - `[["int_se"]]`: SE of this intercept,
#'  - `[["h2"]]`: LDSC regression estimate of (SNP) heritability (also see
#'    [coef_to_liab]),
#'  - `[["h2_se"]]`: SE of this heritability estimate.
#'
#' @importFrom bigassertr assert_one_int
#'
#' @export
#'
snp_ldsc_rg <- function(ld_score, ld_size, z1, z2, sample_size_1, sample_size_2,
                     blocks = 200,
                     intercept = NULL,
                     step1_chisq_max = 30,
                     chi2_thr2 = Inf,
                     ncores = 1, h2_se = FALSE) {

  assert_lengths(z1, z2)

  assert_lengths(z1, ld_score)
  assert_one_int(ld_size)

  M <- length(z1)
  if (length(sample_size_1) == 1) {
    sample_size_1 <- rep(sample_size_1, M)
  } else {
    assert_lengths(sample_size_1, z1)
  }
  if (length(sample_size_2) == 1) {
    sample_size_2 <- rep(sample_size_2, M)
  } else {
    assert_lengths(sample_size_2, z2)
  }


  # First compute heritabilities
  if(h2_se){
    if(is.null(blocks)) stop("If h2_se is desired, please provide blocks.\n")
    h2_blocks <- blocks
  }else{
    h2_blocks <- NULL
  }

  h2_1 <- snp_ldsc(ld_score, ld_size, z1^2,
                   sample_size_1,
                   blocks = h2_blocks,
                   chi2_thr1 = step1_chisq_max,
                   chi2_thr2 = chi2_thr2)
  pred_h2_1 <- h2_1[["int"]] + h2_1[["h2"]]*sample_size_1*ld_score/ld_size


  h2_2 <- snp_ldsc(ld_score, ld_size, z2^2,
                   sample_size_2,
                   blocks = h2_blocks,
                   chi2_thr1 = step1_chisq_max,
                   chi2_thr2 = chi2_thr2)
  pred_h2_2 <- h2_2[["int"]] + h2_2[["h2"]]*sample_size_2*ld_score/ld_size



  step1_index <- which(z1^2 < step1_chisq_max & z2^2 < step1_chisq_max)

  result <- snp_ldsc(ld_score, ld_size, z1*z2,
                     sample_size = sqrt(sample_size_1*sample_size_2),
                     blocks = blocks,
                     intercept= intercept,
                     chi2_thr2 = chi2_thr2,
                     step1_index = step1_index,
                     w0 = pred_h2_1*pred_h2_2,
                     type = "rg")
  gencorr <- result[["h2"]]/sqrt(h2_1[["h2"]]*h2_2[["h2"]])

  res <- c(int = result[["int"]],
          gencov = result[["h2"]],
          t1_int = h2_1[["int"]],
          t2_int = h2_2[["int"]],
          h2_1 = h2_1[["h2"]],
          h2_2 = h2_2[["h2"]],
          gencorr = gencorr)
  if(!is.null(blocks)){
    res <- c(res, int_s = result[["int_se"]],
             gencov_se = result[["h2_se"]])
  }
  if(h2_se){
    res <- c(res,
             t1_int_se = h2_1[["int_se"]],
             t2_int_se  = h2_2[["int_se"]],
             h2_1_se = h2_1[["h2_se"]],
             h2_2_se = h2_2[["h2_se"]])
  }

  return(res)
}

