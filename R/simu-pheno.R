################################################################################

scaled_prod <- function(X, ind, ind.row, ind.col, beta) {
  ms <- big_scale()(X, ind.row, ind.col[ind])
  big_prodVec(X, beta[ind], ind.row, ind.col[ind],
              center = ms$center, scale = ms$scale)
}

#' Simulate phenotypes
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
#'   - "pheno": vector of phenotypes,
#'   - "set": indices of causal variants,
#'   - "effects": effect sizes corresponding to `set`.
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
    # assert_package("rmutil")
    if (!requireNamespace("rmutil", quietly = TRUE))
      stop2("Please install package 'rmutil'.")
    rmutil::rlaplace(M, s = sqrt(h2 / (2 * M)))
  }

  # simulate genetic liability
  y.simu <- big_apply(G, scaled_prod,
                      a.combine = bigparallelr::plus, ind = seq_along(set),
                      ind.row = ind.row, ind.col = set, beta = effects,
                      ncores = ncores)

  # make sure genetic liability has variance equal to heritability
  y.simu <- y.simu / stats::sd(y.simu) * sqrt(h2)
  stopifnot(all.equal(drop(stats::var(y.simu)), h2))

  # add environmental part
  y.simu <- y.simu + stats::rnorm(length(y.simu), sd = sqrt(1 - h2))

  # make binary outcome using liability threshold model
  if (!is.null(K)) y.simu <- as.integer(y.simu > stats::qnorm(1 - K))

  list(pheno = y.simu, set = set, effects = effects)
}

################################################################################
