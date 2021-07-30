################################################################################

#' lassosum2
#'
#' @inheritParams snp_ldpred2_grid
#' @param delta Vector of shrinkage parameters to try (L2-regularization).
#'   Default is `c(0.001, 0.005, 0.02, 0.1, 0.6, 3)`.
#' @param nlambda Number of different lambdas to try (L1-regularization).
#'   Default is `20`.
#' @param lambda.min.ratio Ratio between last and first lambdas to try.
#'   Default is `0.01`.
#' @param dfmax Maximum number of non-zero effects in the model.
#'   Default is `200e3`.
#' @param maxiter Maximum number of iterations before convergence.
#'   Default is `500`.
#' @param tol Tolerance parameter for assessing convergence.
#'   Default is `1e-5`.
#'
#' @return A matrix of effect sizes, one vector (column) for each row in
#'   `attr(<res>, "grid_param")`. Missing values are returned when strong
#'   divergence is detected.
#'
#' @export
#'
snp_lassosum2 <- function(corr, df_beta,
                          delta = signif(seq_log(1e-3, 3, 6), 1),
                          nlambda = 20, lambda.min.ratio = 0.01,
                          dfmax = 200e3, maxiter = 500, tol = 1e-5,
                          ncores = 1) {

  assert_df_with_names(df_beta, c("beta", "beta_se", "n_eff"))
  assert_lengths(rows_along(corr), cols_along(corr), rows_along(df_beta))
  assert_pos(df_beta$beta_se, strict = TRUE)
  assert_pos(delta, strict = TRUE)
  assert_cores(ncores)

  N <- df_beta$n_eff
  scale <- sqrt(N * df_beta$beta_se^2 + df_beta$beta^2)
  beta_hat <- df_beta$beta / scale

  lambda0 <- max(abs(beta_hat))
  seq_lam <- seq_log(lambda0, lambda.min.ratio * lambda0, nlambda + 1)[-1]
  grid_param <- expand.grid(lambda = seq_lam, delta = delta)

  ord <- with(grid_param, order(lambda * (1 + delta)))
  inv_ord <- match(seq_along(ord), ord)

  bigparallelr::register_parallel(ncores)

  res_grid <- foreach(ic = ord, .export = "lassosum2") %dopar% {

    time <- system.time(
      # lassosum2 model
      res <- lassosum2(
        corr     = corr,
        beta_hat = beta_hat,
        lambda   = grid_param$lambda[ic],
        delta    = grid_param$delta[ic],
        dfmax    = dfmax,
        maxiter  = maxiter,
        tol      = tol
      )
    )

    res$time <- time[["elapsed"]]
    res
  }
  res_grid <- res_grid[inv_ord]

  grid_param$num_iter <- sapply(res_grid, function(.) .$num_iter)
  grid_param$time <- sapply(res_grid, function(.) .$time)
  beta_grid <- do.call("cbind", lapply(res_grid, function(.) .$beta_est))
  grid_param$sparsity <- colMeans(beta_grid == 0)

  structure(sweep(beta_grid, 1, scale, '*'), grid_param = grid_param)
}

################################################################################
