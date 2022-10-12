################################################################################

center_vec <- function(x) { x - mean(x) }

################################################################################

#' Simulate phenotypes
#'
#' Simulate phenotypes using a linear model. When a prevalence is given, the
#' liability threshold is used to convert liabilities to a binary outcome.
#' The genetic and environmental liabilities are scaled such that the variance
#' of the genetic liability is exactly equal to the requested heritability, and
#' the variance of the total liability is equal to 1.
#'
#' @inheritParams bigsnpr-package
#' @param h2 Heritability.
#' @param M Number of causal variants.
#' @param ind.possible Indices of possible causal variants.
#' @param effects.dist Distribution of effects.
#'   Either `"gaussian"` (the default) or `"laplace"`.
#' @param K Prevalence. Default is `NULL`, giving a continuous trait.
#' @param alpha Assumes that the average contribution (e.g. heritability)
#'   of a SNP of frequency \eqn{p} is proportional to
#'   \eqn{[2p(1-p)]^{1+\alpha}}. Default is `-1`.
#' @param prob Vector of probability weights for sampling causal indices.
#'   It can have 0s (discarded) and is automatically scaled to sum to 1.
#'   Default is `NULL` (all indices have the same probability).
#'
#' @return A list with 3 elements:
#'   - `$pheno`: vector of phenotypes,
#'   - `$set`: indices of causal variants,
#'   - `$effects`: effect sizes (of scaled genotypes) corresponding to `set`.
#'   - `$allelic_effects`: effect sizes, but on the allele scale (0|1|2).
#' @export
#'
snp_simuPheno <- function(G, h2, M,
                          K = NULL,
                          alpha = -1,
                          ind.row = rows_along(G),
                          ind.possible = cols_along(G),
                          prob = NULL,
                          effects.dist = c("gaussian", "laplace"),
                          ncores = 1) {

  # sample causal variants
  set <- sort(sample(ind.possible, size = M, prob = prob, replace = FALSE))

  # variance of genotypes
  sd <- sqrt(
    big_colstats(G, ind.row = ind.row, ind.col = set, ncores = ncores)$var)

  # sample effect sizes (for causal variants)
  effects <- if (match.arg(effects.dist) == "gaussian") {
    stats::rnorm(M, sd = sd^alpha)
  } else {
    assert_package("rmutil")
    rmutil::rlaplace(M, s = sd^alpha)
  }

  # compute genetic liability
  gen_liab <- big_prodVec(G, effects, ind.row = ind.row, ind.col = set,
                          ncores = ncores)

  # make sure genetic liability has variance equal to heritability
  coeff1 <- sqrt(h2) / stats::sd(gen_liab)
  gen_liab <- center_vec(gen_liab * coeff1)
  stopifnot(all.equal(mean(gen_liab), 0))
  stopifnot(all.equal(stats::var(gen_liab), h2))

  # add environmental part + make sure that total variance is exactly 1
  env_liab <- stats::rnorm(length(gen_liab), sd = sqrt(1 - h2))
  var_env <- stats::var(env_liab)
  cov_env <- stats::cov(gen_liab, env_liab)
  coeff2 <- (sqrt(cov_env^2 + (1 - h2) * var_env) - cov_env) / var_env
  full_liab <- gen_liab + center_vec(env_liab * coeff2)
  stopifnot(all.equal(mean(full_liab), 0))
  stopifnot(all.equal(stats::var(full_liab), 1))

  # possibly make binary outcome using liability threshold model
  pheno <- if (is.null(K)) {
    full_liab
  } else {
    (full_liab > stats::qnorm(K, lower.tail = FALSE)) + 0L
  }

  # return
  list(pheno           = pheno,
       set             = set,
       effects         = effects * coeff1 * sd,
       allelic_effects = effects * coeff1)
}

################################################################################
