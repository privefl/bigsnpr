#' Title
#'
#' @param x
#' @param covar
#' @param ncores
#' @param tol
#' @param maxiter
#'
#' @return
#' @export
#'
#' @examples
GWAS <- function(x, covar = NULL, ncores = 1, tol = 1e-8, maxiter = 100) {
  check_x(x, check.y = TRUE)

  X <- x$genotypes
  X.desc <- describe(X)
  y <- (x$fam$pheno + 1) / 2

  if (is.null(covar)) {
    n <- nrow(X)
    covar <- cbind(rep(0, n), rep(1, n))
  } else {
    covar <- cbind(0, 1, covar)
  }

  mod0 = glm(y ~ covar - 1, family = binomial)
  p0 = mod0$fitted
  w0 = p0 * (1 - p0)
  z0 = log(p0 / (1 - p0)) + (y - p0) / w0
  rm(mod0, p0)

  GWAS.part <- function(lims) {
    X.part <- sub.big.matrix(X.desc,
                             firstCol = lims[1],
                             lastCol = lims[2],
                             backingpath = x$backingpath)

    res <- wcrossprod(X.part@address, covar, y, z0, w0,
                      tol, maxiter)

    indNoConv <- which(!res$conv)
    if ((l <- length(indNoConv)) > 0)
      printf(paste("For %d SNPs, IRLS has not converged",
                   "using glm for those instead.\n", sep = "; "), l)
    for (j in indNoConv) {
      mod <- glm(y ~ X.part[, j] + covar - 1, family = binomial)
      coeffs <- summary(mod)$coefficients
      res$betas[j] <- coeffs[1]
      res$std[j] <- coeffs[2]
    }

    rbind(res$betas, res$std)
  }

  range.parts <- CutBySize(ncol(X), nb = ncores)

  obj <- foreach::foreach(i = 1:nrow(range.parts),
                          .combine = "cbind",
                          .noexport = c('x', 'X'))
  expr_fun <- function(i) GWAS.part(range.parts[i, ])
  foreach2(obj, expr_fun, ncores)
}
