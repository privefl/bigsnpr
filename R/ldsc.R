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

  assert_pos(chi2)
  assert_lengths(chi2, ld_div_size)
  if (length(sample_size) != 1) assert_lengths(sample_size, chi2)
  if (length(blocks) == 1) blocks <- sort(rep_len(seq_len(blocks), length(chi2)))

  #### step 1 ####

  ind_sub1 <- which(chi2 < chi2_thr1)
  x1 <- (ld_div_size * sample_size)[ind_sub1]
  y1 <- chi2[ind_sub1]
  ind_blocks <- split(seq_along(ind_sub1), blocks[ind_sub1])
  n_blocks <- length(ind_blocks)

  pred0 <- y1
  for (i in 1:100) {
    w1 <- WEIGHTS(pred0, x1)
    pred <- stats::lm(y1 ~ x1, weights = w1)$fitted.values
    if (max(abs(pred - pred0)) < 1e-6) break
    pred0 <- pred
  }
  xw1 <- WEIGHT(w1, cbind(x1, 1))
  yw1 <- WEIGHT(w1, y1)

  xtx_block_values <- xty_block_values <- list()
  for (i in seq_along(ind_blocks)) {
    ind_block_i <- ind_blocks[[i]]
    X <- xw1[ind_block_i, , drop = FALSE]
    xtx_block_values[[i]] <- crossprod(X)
    xty_block_values[[i]] <- crossprod(X, yw1[ind_block_i])
  }
  xtx <- Reduce('+', xtx_block_values)
  xty <- Reduce('+', xty_block_values)

  delete_values <- sapply(seq_len(n_blocks), function(i) {
    solve(xtx - xtx_block_values[[i]], xty - xty_block_values[[i]])[2]
  })
  pseudovalues <- n_blocks * solve(xtx, xty)[2] - (n_blocks - 1) * delete_values

  step1_int    <- mean(pseudovalues)
  step1_int_se <- sd(pseudovalues / n_blocks)

  #### step 2 ####

  ind_sub2 <- which(chi2 < chi2_thr2)
  x <- (ld_div_size * sample_size)[ind_sub2]
  y <- chi2[ind_sub2]
  yp <- y - step1_int
  ind_blocks <- split(seq_along(ind_sub2), blocks[ind_sub2])
  n_blocks <- length(ind_blocks)

  pred0 <- y
  for (i in 1:100) {
    w2 <- WEIGHTS(pred0, x)
    pred <- step1_int + stats::lm(yp ~ x + 0, weights = w2)$fitted.values
    if (max(abs(pred - pred0)) < 1e-6) break
    pred0 <- pred
  }
  xw2 <- WEIGHT(w2, x)
  yw2 <- WEIGHT(w2, yp)

  xtx_block_values <- xty_block_values <- list()
  for (i in seq_along(ind_blocks)) {
    ind_block_i <- ind_blocks[[i]]
    X <- xw2[ind_block_i]
    xtx_block_values[[i]] <- crossprod(X)
    xty_block_values[[i]] <- crossprod(X, yw2[ind_block_i])
  }
  xtx <- Reduce('+', xtx_block_values)
  xty <- Reduce('+', xty_block_values)

  delete_values <- sapply(seq_len(n_blocks), function(i) {
    solve(xtx - xtx_block_values[[i]], xty - xty_block_values[[i]])
  })
  pseudovalues <- n_blocks * c(solve(xtx, xty)) - (n_blocks - 1) * delete_values

  c(int    = step1_int,
    int_se = step1_int_se,
    h2     = mean(pseudovalues),
    h2_se  = sd(pseudovalues / n_blocks))
}

################################################################################
