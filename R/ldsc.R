################################################################################

# heteroscedasticity and overcounting weights
WEIGHTS <- function(pred, w_ld) {
  1 / (pred^2 * w_ld)
}

crossprod2 <- function(x, y) drop(base::crossprod(x, y))

# equivalent to stats::lm.wfit(cbind(1, x), y, w)
wlm <- function(x, y, w) {
  wx <- w * x
  W   <- sum(w)
  WX  <- sum(wx)
  WY  <- crossprod2(w,  y)
  WXX <- crossprod2(wx, x)
  WXY <- crossprod2(wx, y)
  alpha <- (WXX * WY - WX * WXY) / (W * WXX - WX^2)
  beta  <- (WXY * W  - WX * WY)  / (W * WXX - WX^2)
  list(intercept = alpha, slope = beta, pred = x * beta + alpha)
}

# equivalent to stats::lm.wfit(as.matrix(x), y, w)
wlm_no_int <- function(x, y, w) {
  wx <- w * x
  WXX <- crossprod2(wx, x)
  WXY <- crossprod2(wx, y)
  beta  <- WXY / WXX
  list(slope = beta, pred = x * beta)
}

################################################################################

#' LD score regression
#'
#' Implementation of LDSC regression in R, originally based on
#' https://helda.helsinki.fi/bitstream/handle/10138/273501/MT_ldsc_hautakangas.pdf
#'
#' @param ld_score Vector of LD scores.
#' @param ld_size Number of variants used to compute `ld_score`.
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
#' @importFrom bigassertr assert_one_int
#'
#' @export
#'
snp_ldsc <- function(ld_score, ld_size, chi2, sample_size, blocks = 200,
                     chi2_thr1 = 30, chi2_thr2 = Inf) {

  assert_pos(chi2, strict = FALSE)
  assert_lengths(chi2, ld_score)
  assert_one_int(ld_size)
  if (length(sample_size) == 1) {
    sample_size <- rep(sample_size, length(chi2))
  } else {
    assert_lengths(sample_size, chi2)
  }

  if (is.null(blocks)) {

    #### step 1 ####

    ind_sub1 <- which(chi2 < chi2_thr1)
    w_ld <- pmax(ld_score[ind_sub1], 1)
    x1 <- (ld_score / ld_size * sample_size)[ind_sub1]
    y1 <- chi2[ind_sub1]

    pred0 <- y1
    for (i in 1:100) {
      pred <- wlm(x1, y1, w = WEIGHTS(pred0, w_ld))$pred
      if (max(abs(pred - pred0)) < 1e-6) break
      pred0 <- pred
    }
    step1_int <- wlm(x1, y1, w = WEIGHTS(pred0, w_ld))$intercept

    #### step 2 ####

    ind_sub2 <- which(chi2 < chi2_thr2)
    w_ld <- pmax(ld_score[ind_sub2], 1)
    x <- (ld_score / ld_size * sample_size)[ind_sub2]
    y <- chi2[ind_sub2]
    yp <- y - step1_int

    pred0 <- y
    for (i in 1:100) {
      pred <- step1_int + wlm_no_int(x, yp, w = WEIGHTS(pred0, w_ld))$pred
      if (max(abs(pred - pred0)) < 1e-6) break
      pred0 <- pred
    }
    step2_h2 <- wlm_no_int(x, yp, w = WEIGHTS(pred0, w_ld))$slope

    c(int = step1_int, h2 = step2_h2)

  } else {

    est <- snp_ldsc(ld_score, ld_size, chi2, sample_size,
                    NULL, chi2_thr1, chi2_thr2)

    if (length(blocks) == 1)
      blocks <- sort(rep_len(seq_len(blocks), length(chi2)))

    ind_blocks <- split(seq_along(blocks), blocks)
    n_blocks <- length(ind_blocks)

    delete_values <- sapply(1:n_blocks, function(k) {
      ind_rm <- ind_blocks[[k]]
      snp_ldsc(ld_score[-ind_rm], ld_size, chi2[-ind_rm], sample_size[-ind_rm],
               NULL, chi2_thr1, chi2_thr2)
    })
    pseudovalues <- n_blocks * est - (n_blocks - 1) * delete_values

    c(int    = mean(pseudovalues[1, ]),
      int_se = sd(pseudovalues[1, ]) / sqrt(n_blocks),
      h2     = mean(pseudovalues[2, ]),
      h2_se  = sd(pseudovalues[2, ]) / sqrt(n_blocks))

  }
}

################################################################################
