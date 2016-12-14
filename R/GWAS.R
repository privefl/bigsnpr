#' Title
#'
#' @param x
#' @param covar
#' @param ncores
#' @param tol
#' @param maxiter
#'
#' @import foreach
#'
#' @return
#' @export
#'
#' @examples
snp_logitGWAS <- function(x, ind.train = seq(nrow(X)),
                       covar = NULL, ncores = 1,
                       tol = 1e-8, maxiter = 100) {
  check_x(x)

  X <- x$genotypes
  X.desc <- describe(X)
  y.train <- transform_levels(x$fam$affection[ind.train])
  n <- length(ind.train)

  if (is.null(covar)) {
    covar <- cbind(rep(0, n), rep(1, n))
  } else {
    covar <- cbind(0, 1, covar)
  }
  stopifnot(n == nrow(covar))

  # no intercept because already in covar
  mod0 <- stats::glm(y.train ~ covar - 1, family = binomial)
  p0 <- mod0$fitted
  w0 <- p0 * (1 - p0)
  z0 <- log(p0 / (1 - p0)) + (y.train - p0) / w0
  rm(mod0, p0)

  PATH <- x$backingpath
  range.parts <- CutBySize(ncol(X), nb = ncores)

  if (is.seq <- (ncores == 1)) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  }
  res.all <- foreach(ic = seq_len(ncores), .combine = 'cbind') %dopar% {
    lims <- range.parts[ic, ]

    X.part <- sub.big.matrix(X.desc,
                             firstCol = lims[1],
                             lastCol = lims[2],
                             backingpath = PATH)

    # https://www.r-bloggers.com/too-much-parallelism-is-as-bad/
    multi <- (!is.seq) && detect_MRO()
    if (multi) nthreads.save <- RevoUtilsMath::setMKLthreads(1)
    res <- wcrossprod(X.part@address, covar, y.train, z0, w0,
                      ind.train, tol, maxiter)
    if (multi) RevoUtilsMath::setMKLthreads(nthreads.save)

    indNoConv <- which(!res$conv)
    if ((l <- length(indNoConv)) > 0) {
      printf(paste("For %d SNPs, IRLS has not converged",
                   "using glm for those instead.\n", sep = "; "), l)

      for (j in indNoConv) {
        mod <- glm(y.train ~ X.part[ind.train, j] + covar - 1, family = binomial)
        coeffs <- summary(mod)$coefficients
        res$betas[j] <- coeffs[1]
        res$std[j] <- coeffs[2]
      }
    }

    rbind(res$betas, res$std)
  }
  if (!is.seq) parallel::stopCluster(cl)

  res.all
}
