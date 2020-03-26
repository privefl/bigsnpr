################################################################################

WEIGHTS <- function(w_ld, pred) {
  het_w <- 1 / (2 * pred^2)       # heteroscedasticity weights
  oc_w  <- 1 / pmax(w_ld, 1)      # overcounting weights
  het_w * oc_w                    # merged weights
}

WEIGHT <- function(w, x) {
  w1 <- sqrt(w)
  w_norm <- w1 / sum(w1)
  x * w_norm
}

################################################################################

#' LD score regression
#'
#' Implementation of LDSC regression in R, originally based on
#' https://helda.helsinki.fi/bitstream/handle/10138/273501/MT_ldsc_hautakangas.pdf
#'
#' @param ld_div_size Vector of LD scores, divided by the number of variants
#'   used to compute these LD scores.
#' @param chi2 Vector of chi-squared statistics.
#' @param sample_size Sample size of GWAS corresponding to chi-squared statistics.
#'   Possibly a vector, or just a single value.
#' @param chi2_thr1 Threshold on `chi2` in step 1. Default is `30`.
#' @param chi2_thr2 Threshold on `chi2` in step 2. Default is `Inf` (none).
#' @param blocks Either a simgle number specifying the number of blocks,
#'   or a vector of integers specifying the block number of each `chi2` value.
#'   Default is `200`, dividing into 200 blocks of approximately equal size.
#'
#' @return Vector of 4 values:
#'  1. LDSC regression intercept,
#'  2. SE of this intercept,
#'  3. LDSC regression estimate of heritability,
#'  4. SE of this heritability estimate.
#'
#' @export
#'
snp_ldsc <- function(ld_div_size, chi2, sample_size, blocks = 200,
                     chi2_thr1 = 30, chi2_thr2 = Inf) {

  assert_pos(chi2, strict = FALSE)
  assert_lengths(chi2, ld_div_size)
  if (length(sample_size) == 1) {
    sample_size <- rep(sample_size, length(chi2))
  } else {
    assert_lengths(sample_size, chi2)
  }

  if (is.null(blocks)) {

    #### step 1 ####

    ind_sub1 <- which(chi2 < chi2_thr1)
    x1 <- (ld_div_size * sample_size)[ind_sub1]
    y1 <- chi2[ind_sub1]

    pred0 <- y1
    for (i in 1:100) {
      w1 <- WEIGHTS(pred0, x1)
      pred <- stats::lm(y1 ~ x1, weights = w1)$fitted.values
      if (max(abs(pred - pred0)) < 1e-6) break
      pred0 <- pred
    }
    xw1 <- WEIGHT(w1, cbind(x1, 1))
    yw1 <- WEIGHT(w1, y1)

    xtx <- crossprod(xw1)
    xty <- crossprod(xw1, yw1)
    step1_int <- solve(xtx, xty)[2]


    #### step 2 ####

    ind_sub2 <- which(chi2 < chi2_thr2)
    x <- (ld_div_size * sample_size)[ind_sub2]
    y <- chi2[ind_sub2]
    yp <- y - step1_int

    pred0 <- y
    for (i in 1:100) {
      w2 <- WEIGHTS(pred0, x)
      pred <- step1_int + stats::lm(yp ~ x + 0, weights = w2)$fitted.values
      if (max(abs(pred - pred0)) < 1e-6) break
      pred0 <- pred
    }
    xw2 <- WEIGHT(w2, x)
    yw2 <- WEIGHT(w2, yp)

    xtx <- crossprod(xw2)
    xty <- crossprod(xw2, yw2)
    step2_h2 <- solve(xtx, xty)[1]

    c(int = step1_int, h2 = step2_h2)

  } else {

    est <- snp_ldsc(ld_div_size, chi2, sample_size, NULL,
                    chi2_thr1, chi2_thr2)

    if (length(blocks) == 1)
      blocks <- sort(rep_len(seq_len(blocks), length(chi2)))

    ind_blocks <- split(seq_along(blocks), blocks)
    n_blocks <- length(ind_blocks)

    delete_values <- sapply(1:n_blocks, function(k) {
      ind <- unlist(ind_blocks[-k])
      snp_ldsc(ld_div_size[ind], chi2[ind], sample_size[ind], NULL,
               chi2_thr1, chi2_thr2)
    })
    pseudovalues <- n_blocks * est - (n_blocks - 1) * delete_values

    c(int    = mean(pseudovalues[1, ]),
      int_se = sd(pseudovalues[1, ]) / sqrt(n_blocks),
      h2     = mean(pseudovalues[2, ]),
      h2_se  = sd(pseudovalues[2, ]) / sqrt(n_blocks))

  }
}

################################################################################
