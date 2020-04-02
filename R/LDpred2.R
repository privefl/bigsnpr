################################################################################

#' LDpred2
#'
#' @inheritParams snp_ldsc2
#' @param h2 Heritability estimate.
#'   Default is estimated using constrained LD score regression.
#'
#' @return `snp_ldpred2_inf`: A vector of effects, assuming an infinitesimal model.
#'
#' @export
#'
#' @rdname LDpred2
#'
snp_ldpred2_inf <- function(corr, df_beta, h2 = NULL) {

  assert_lengths(rows_along(corr), cols_along(corr), rows_along(df_beta))
  stopifnot(all.equal(colnames(df_beta), c("beta", "beta_se", "n_eff")))

  if (is.null(h2))
    h2 <- snp_ldsc2(corr, df_beta, intercept = 1, blocks = NULL)[["h2"]]

  m <- ncol(corr)
  N <- df_beta$n_eff

  sd <- df_beta$beta_se * sqrt(N)
  beta_hat <- df_beta$beta / sd
  beta_inf <- as.vector(Matrix::solve(
    corr + Matrix::Diagonal(m, m / (h2 * N)), beta_hat))

  beta_inf * sd
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
                             burn_in = 100,
                             num_iter = 200,
                             ncores = 1) {

  assert_lengths(rows_along(corr), cols_along(corr), rows_along(df_beta))
  stopifnot(all.equal(colnames(df_beta), c("beta", "beta_se", "n_eff")))
  stopifnot(all.equal(colnames(grid_param), c("p", "h2", "sparse")))

  N <- df_beta$n_eff
  sd <- df_beta$beta_se * sqrt(N)
  beta_hat <- df_beta$beta / sd

  # compute one infinitesimal model, just for initialization
  m <- ncol(corr)
  beta_inf <- as.vector(Matrix::solve(
    corr + Matrix::Diagonal(m, m / (median(grid_param$h2) * N)), beta_hat))

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

  sweep(beta_gibbs, 1, sd, '*')
}

################################################################################

#' @param vec_p_init Vector of initial values for p. Default is `0.1`.
#'   They can be run in parallel by changing `ncores`.
#' @param h2_init Heritability estimate for initialization.
#'   Default is estimated using constrained LD score regression.
#' @param verbose Whether to print "p // h2" estimates at each iteration.
#'
#' @return `snp_ldpred2_auto`: A list of `length(vec_p_init)` list(s) with
#'   - `$beta_est`: vector of effect sizes
#'   - `$p_est`: estimate of p, the proportion of causal variants
#'   - `$h2_est`: estimate of the (SNP) heritability
#'   - `$path_p_est`: full path of p estimates (including burn-in);
#'     useful to check convergence of the iterative algorithm
#'   - `$path_h2_est`: full path of h2 estimates (including burn-in);
#'     useful to check convergence of the iterative algorithm
#'
#' @export
#'
#' @rdname LDpred2
#'
snp_ldpred2_auto <- function(corr, df_beta,
                             vec_p_init = 0.1,
                             h2_init = NULL,
                             burn_in = 2000,
                             num_iter = 500,
                             verbose = FALSE,
                             ncores = 1) {

  assert_lengths(rows_along(corr), cols_along(corr), rows_along(df_beta))
  stopifnot(all.equal(colnames(df_beta), c("beta", "beta_se", "n_eff")))
  stopifnot(!(verbose && (ncores > 1)))

  N <- df_beta$n_eff
  sd <- df_beta$beta_se * sqrt(N)
  beta_hat <- df_beta$beta / sd

  if (is.null(h2_init))
    h2_init <- snp_ldsc2(corr, df_beta, intercept = 1, blocks = NULL)[["h2"]]

  # compute one infinitesimal model, just for initialization
  m <- ncol(corr)
  beta_inf <- as.vector(Matrix::solve(
    corr + Matrix::Diagonal(m, m / (h2_init * N)), beta_hat))

  all_ldpred_auto <- ldpred2_gibbs_auto(
    corr      = corr,
    beta_hat  = beta_hat,
    beta_init = beta_inf,
    order     = order(beta_inf^2, decreasing = TRUE) - 1L,
    n_vec     = N,
    p_init    = vec_p_init,
    burn_in   = burn_in,
    num_iter  = num_iter,
    verbose   = verbose,
    ncores    = ncores
  )

  for (k in seq_along(all_ldpred_auto)) {
    all_ldpred_auto[[k]]$beta_est <- drop(all_ldpred_auto[[k]]$beta_est) * sd
  }
  all_ldpred_auto
}

################################################################################
