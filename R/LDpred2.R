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

  assert_df_with_names(df_beta, c("beta", "beta_se", "n_eff"))
  assert_lengths(rows_along(corr), cols_along(corr), rows_along(df_beta))
  assert_pos(df_beta$beta_se, strict = TRUE)

  if (is.null(h2))
    h2 <- snp_ldsc2(corr, df_beta, intercept = 1, blocks = NULL)[["h2"]]
  assert_pos(h2, strict = TRUE)

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
#'   for each row of `grid_param`. When the algorithm diverges, 0s are returned.
#'
#' @export
#'
#' @rdname LDpred2
#'
snp_ldpred2_grid <- function(corr, df_beta, grid_param,
                             burn_in = 100,
                             num_iter = 200,
                             ncores = 1,
                             tmpdir = tempdir()) {

  assert_df_with_names(df_beta, c("beta", "beta_se", "n_eff"))
  assert_df_with_names(grid_param, c("p", "h2", "sparse"))
  assert_lengths(rows_along(corr), cols_along(corr), rows_along(df_beta))
  assert_pos(df_beta$beta_se, strict = TRUE)
  assert_cores(ncores)

  N <- df_beta$n_eff
  sd <- df_beta$beta_se * sqrt(N)
  beta_hat <- df_beta$beta / sd

  # compute one infinitesimal model, just for initialization
  m <- ncol(corr)
  assert_pos(grid_param$h2, strict = TRUE)
  h2_init <- stats::median(grid_param$h2)
  beta_inf <- as.vector(Matrix::solve(
    corr + Matrix::Diagonal(m, m / (h2_init * N)), beta_hat))

  tmp <- tempfile(tmpdir = tmpdir)

  beta_gibbs <- ldpred2_gibbs(
    corr      = bigsparser::as_SFBM(as(corr, "dgCMatrix"), backingfile = tmp),
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

  file.remove(paste0(tmp, ".sbk"))

  sweep(beta_gibbs, 1, sd, '*')
}

################################################################################

#' @param p_init Initial value for p. Default is `0.1`.
#' @param h2_init Heritability estimate for initialization.
#'   Default is estimated using constrained LD score regression.
#' @param verbose Whether to print "p // h2" estimates at each iteration.
#'
#' @return `snp_ldpred2_auto`: A list with
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
                             p_init = 0.1,
                             h2_init = NULL,
                             burn_in = 1000,
                             num_iter = 1000,
                             verbose = FALSE) {

  assert_df_with_names(df_beta, c("beta", "beta_se", "n_eff"))
  assert_lengths(rows_along(corr), cols_along(corr), rows_along(df_beta))
  assert_pos(df_beta$beta_se, strict = TRUE)

  N <- df_beta$n_eff
  sd <- df_beta$beta_se * sqrt(N)
  beta_hat <- df_beta$beta / sd

  if (is.null(h2_init))
    h2_init <- snp_ldsc2(corr, df_beta, intercept = 1, blocks = NULL)[["h2"]]
  assert_pos(h2_init, strict = TRUE)

  # compute one infinitesimal model, just for initialization
  m <- ncol(corr)
  beta_inf <- as.vector(Matrix::solve(
    corr + Matrix::Diagonal(m, m / (h2_init * N)), beta_hat))

  ldpred_auto <- ldpred2_gibbs_auto(
    corr      = corr,
    beta_hat  = beta_hat,
    beta_init = beta_inf,
    order     = order(beta_inf^2, decreasing = TRUE) - 1L,
    n_vec     = N,
    p_init    = p_init,
    burn_in   = burn_in,
    num_iter  = num_iter,
    verbose   = verbose
  )
  ldpred_auto$beta_est <- drop(ldpred_auto$beta_est) * sd

  ldpred_auto
}

################################################################################
