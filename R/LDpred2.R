################################################################################

#' LDpred2
#'
#' LDpred2. Tutorial at \url{https://bit.ly/ldpred2}.
#'
#' @inheritParams snp_ldsc2
#' @param corr Sparse correlation matrix as an [SFBM][SFBM-class]. If `corr` is
#'   a dgSMatrix, you can use `bigsparser::as_SFBM(as(corr, "dgCMatrix"))`.
#' @param h2 Heritability estimate.
#'
#' @return `snp_ldpred2_inf`: A vector of effects, assuming an infinitesimal model.
#'
#' @export
#'
#' @rdname LDpred2
#'
snp_ldpred2_inf <- function(corr, df_beta, h2) {

  assert_df_with_names(df_beta, c("beta", "beta_se", "n_eff"))
  assert_lengths(rows_along(corr), cols_along(corr), rows_along(df_beta))
  assert_pos(df_beta$beta_se, strict = TRUE)
  assert_pos(h2, strict = TRUE)

  N <- df_beta$n_eff
  scale <- df_beta$beta_se * sqrt(N)
  beta_hat <- df_beta$beta / scale

  beta_inf <- bigsparser::sp_solve_sym(
    corr, beta_hat, add_to_diag = ncol(corr) / (h2 * N))

  beta_inf * scale
}

################################################################################

#' @param grid_param A data frame with 3 columns as a grid of hyper-parameters:
#'   - `$p`: proportion of causal variants
#'   - `$h2`: heritability (captured by the variants used)
#'   - `$sparse`: boolean, whether a sparse model is sought
#'   They can be run in parallel by changing `ncores`.
#' @param burn_in Number of burn-in iterations.
#' @param num_iter Number of iterations after burn-in.
#' @inheritParams bigsnpr-package
#'
#' @return `snp_ldpred2_grid`: A matrix of effect sizes, one vector (column)
#'   for each row of `grid_param`.
#'
#' @export
#'
#' @rdname LDpred2
#'
snp_ldpred2_grid <- function(corr, df_beta, grid_param,
                             burn_in = 50,
                             num_iter = 100,
                             ncores = 1) {

  assert_df_with_names(df_beta, c("beta", "beta_se", "n_eff"))
  assert_df_with_names(grid_param, c("p", "h2", "sparse"))
  assert_lengths(rows_along(corr), cols_along(corr), rows_along(df_beta))
  assert_pos(df_beta$beta_se, strict = TRUE)
  assert_pos(grid_param$h2, strict = TRUE)
  assert_cores(ncores)

  N <- df_beta$n_eff
  scale <- df_beta$beta_se * sqrt(N)
  beta_hat <- df_beta$beta / scale

  # compute one infinitesimal model, just for initialization
  med_h2 <- stats::median(grid_param$h2)
  beta_inf <- bigsparser::sp_solve_sym(
    corr, beta_hat, add_to_diag = ncol(corr) / (med_h2 * N))

  beta_gibbs <- ldpred2_gibbs(
    corr      = corr,
    beta_hat  = beta_hat,
    beta_init = beta_inf,
    order     = order(beta_inf^2, decreasing = TRUE) - 1L,
    n_vec     = N,
    h2        = grid_param$h2,
    p         = grid_param$p,
    sparse    = grid_param$sparse,
    burn_in   = burn_in,
    num_iter  = num_iter,
    ncores    = ncores
  )

  sweep(beta_gibbs, 1, scale, '*')
}

################################################################################

#' @param vec_p_init Vector of initial values for p. Default is `0.1`.
#' @param h2_init Heritability estimate for initialization.
#' @param sparse In LDpred2-auto, whether to also report a sparse solution by
#'   running LDpred2-grid with the estimates of p and h2 from LDpred2-auto, and
#'   sparsity enabled. Default is `FALSE`.
#' @param verbose Whether to print "p // h2" estimates at each iteration.
#'
#' @return `snp_ldpred2_auto`: A list (over `vec_p_init`) of lists with
#'   - `$beta_est`: vector of effect sizes
#'   - `$p_est`: estimate of p, the proportion of causal variants
#'   - `$h2_est`: estimate of the (SNP) heritability (also see [coef_to_liab])
#'   - `$path_p_est`: full path of p estimates (including burn-in);
#'     useful to check convergence of the iterative algorithm
#'   - `$path_h2_est`: full path of h2 estimates (including burn-in);
#'     useful to check convergence of the iterative algorithm
#'
#' @import foreach
#'
#' @export
#'
#' @rdname LDpred2
#'
snp_ldpred2_auto <- function(corr, df_beta, h2_init,
                             vec_p_init = 0.1,
                             burn_in = 1000,
                             num_iter = 500,
                             sparse = FALSE,
                             verbose = FALSE,
                             ncores = 1) {

  assert_df_with_names(df_beta, c("beta", "beta_se", "n_eff"))
  assert_lengths(rows_along(corr), cols_along(corr), rows_along(df_beta))
  assert_pos(df_beta$beta_se, strict = TRUE)
  assert_pos(h2_init, strict = TRUE)

  N <- df_beta$n_eff
  scale <- df_beta$beta_se * sqrt(N)
  beta_hat <- df_beta$beta / scale

  # compute one infinitesimal model, just for initialization
  beta_inf <- bigsparser::sp_solve_sym(
    corr, beta_hat, add_to_diag = ncol(corr) / (h2_init * N))

  bigparallelr::register_parallel(ncores)

  foreach(p_init = vec_p_init) %dopar% {

    ldpred_auto <- ldpred2_gibbs_auto(
      corr      = corr,
      beta_hat  = beta_hat,
      beta_init = beta_inf,
      order     = order(beta_inf^2, decreasing = TRUE) - 1L,
      n_vec     = N,
      p_init    = p_init,
      h2_init   = h2_init,
      burn_in   = burn_in,
      num_iter  = num_iter,
      verbose   = verbose
    )
    ldpred_auto$beta_est <- drop(ldpred_auto$beta_est) * scale

    if (sparse) {
      beta_gibbs <- ldpred2_gibbs(
        corr      = corr,
        beta_hat  = beta_hat,
        beta_init = beta_inf,
        order     = order(beta_inf^2, decreasing = TRUE) - 1L,
        n_vec     = N,
        h2        = ldpred_auto$h2_est,
        p         = ldpred_auto$p_est,
        sparse    = TRUE,
        burn_in   = 50,
        num_iter  = 100,
        ncores    = 1
      )
      ldpred_auto$beta_est_sparse <- drop(beta_gibbs) * scale
    }

    ldpred_auto
  }
}

################################################################################
