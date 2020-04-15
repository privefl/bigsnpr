################################################################################

scaled_prod <- function(X, ind, ind.row, ind.col, beta) {
  ms <- big_scale()(X, ind.row, ind.col[ind])
  big_prodVec(X, beta[ind], ind.row, ind.col[ind],
              center = ms$center, scale = ms$scale)
}

#' Simulate phenotypes
#'
#' Simulate phenotypes using a linear model. When a prevalence is given, the
#' liability threshold is used to convert liabilities to a binary outcome.
#' The genetic and environmental liabilities are scaled such that the variance
#' of the genetic liability is equality the requested heritability, and the
#' variance of the total liability is 1.
#'
#' @inheritParams bigsnpr-package
#' @param h2 Heritability.
#' @param M Number of causal variants.
#' @param ind.possible Indices of possible causal variants.
#' @param effects.dist Distribution of effects.
#'   Either `"gaussian"` (the default) or `"laplace"`.
#' @param K Prevalence. Default is `NULL`, giving a continuous trait.
#'
#' @return A list with 3 elements:
#'   - `$pheno`: vector of phenotypes,
#'   - `$set`: indices of causal variants,
#'   - `$effects`: effect sizes corresponding to `set`.
#' @export
#'
snp_simuPheno <- function(G, h2, M, K = NULL,
                          ind.row = rows_along(G),
                          ind.possible = cols_along(G),
                          effects.dist = c("gaussian", "laplace"),
                          ncores = 1) {

  set <- sample(ind.possible, size = M)
  effects <- if (match.arg(effects.dist) == "gaussian") {
    stats::rnorm(M, sd = sqrt(h2 / M))
  } else {
    assert_package("rmutil")
    rmutil::rlaplace(M, s = sqrt(h2 / (2 * M)))
  }

  # simulate genetic liability
  gen_liab <- big_apply(G, scaled_prod,
                        a.combine = bigparallelr::plus, ind = seq_along(set),
                        ind.row = ind.row, ind.col = set, beta = effects,
                        ncores = ncores)

  # make sure genetic liability has variance equal to heritability
  coeff1 <- sqrt(h2) / stats::sd(gen_liab)
  gen_liab <- gen_liab * coeff1
  stopifnot(all.equal(stats::var(gen_liab), h2))

  # add environmental part + make sure that total variance is 1
  env_liab <- stats::rnorm(length(gen_liab), sd = sqrt(1 - h2))
  var_env <- stats::var(env_liab)
  cov_env <- stats::cov(gen_liab, env_liab)
  coeff2 <- (sqrt(cov_env^2 + (1 - h2) * var_env) - cov_env) / var_env
  full_liab <- gen_liab + env_liab * coeff2
  stopifnot(all.equal(stats::var(full_liab), 1))

  # make binary outcome using liability threshold model
  pheno <- if (is.null(K)) full_liab else (full_liab > stats::qnorm(1 - K)) + 0L

  list(pheno = pheno, set = set, effects = effects * coeff1)
}

################################################################################
