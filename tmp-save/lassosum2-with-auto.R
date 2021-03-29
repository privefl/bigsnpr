################################################################################

#' lassosum2
#'
#' @inheritParams snp_ldpred2_grid
#' @param s Vector of shrinkage parameters to try. Default is `1:10 / 10`.
#'   Try to avoid using `0` as it often leads to effects exploding.
#' @param nlambda Number of different lambdas to try. Default is `20`.
#' @param lambda.min.ratio Ratio between last and first lambdas to try.
#'   Default is `0.01`.
#' @param dfmax Maximum number of non-zero effects in the model.
#'   Default is `200e3`. This is not used when `s = 1`.
#' @param maxiter Maximum number of iterations before convergence.
#'   Default is `500`.
#' @param tol Tolerance parameter for assessing convergence. Default is `1e-5`.
#'
#' @return A matrix of effect sizes, one vector (column) for each row in
#'   `attr(<res>, "grid_param")`, where you can get extra information.
#'
#' @export
#'
snp_lassosum2 <- function(corr, df_beta, s = 1:10 / 10,
                          nlambda = 20, lambda.min.ratio = 0.01,
                          dfmax = 200e3, maxiter = 500, tol = 1e-5,
                          ncores = 1) {

  assert_df_with_names(df_beta, c("beta", "beta_se", "n_eff"))
  assert_lengths(rows_along(corr), cols_along(corr), rows_along(df_beta))
  assert_pos(df_beta$beta_se, strict = TRUE)
  assert_cores(ncores)
  assert_package("fdrtool")

  N <- df_beta$n_eff
  scale <- sqrt(N * df_beta$beta_se^2 + df_beta$beta^2)
  beta_hat <- df_beta$beta / scale

  chi2 <- (df_beta$beta / df_beta$beta_se)^2
  pval <- stats::pchisq(chi2, df = 1, lower.tail = FALSE)
  utils::capture.output(
    fdr <- fdrtool::fdrtool(pval, statistic = "pvalue", plot = FALSE))
  beta_hat_shrunk <- beta_hat * (1 - fdr$lfdr)

  lambda0 <- max(abs(beta_hat))
  seq_lam <- head(seq_log(lambda.min.ratio * lambda0, lambda0, nlambda + 1), -1)
  grid_param <- expand.grid(lambda = seq_lam, s = sort(s))
  res_grid <- FBM(nrow(grid_param), 3) # auto_score / num_iter / time

  bigparallelr::register_parallel(ncores)

  beta_grid <- foreach(ic = rows_along(grid_param), .export = "lassosum2") %dopar% {

    # lassosum2 model
    time <- system.time(
      res <- lassosum2(
        corr     = corr,
        beta_hat = beta_hat,
        lambda   = grid_param$lambda[ic],
        s        = grid_param$s[ic],
        dfmax    = dfmax,
        maxiter  = maxiter,
        tol      = tol
      )
    )

    beta <- res$beta_est
    bRb <- crossprod(beta, bigsparser::sp_prodVec(corr, beta))

    res_grid[ic, 1] <- crossprod(beta, beta_hat_shrunk) / sqrt(bRb)
    res_grid[ic, 2] <- res$num_iter
    res_grid[ic, 3] <- time[[3]]

    beta * scale
  }

  beta_grid <- do.call("cbind", beta_grid)
  grid_param$sparsity   <- colMeans(beta_grid == 0)
  grid_param$auto_score <- res_grid[, 1]
  grid_param$num_iter   <- res_grid[, 2]
  grid_param$time       <- res_grid[, 3]

  structure(beta_grid, grid_param = grid_param)
}

################################################################################
