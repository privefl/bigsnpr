#' @title LD pruning for a "bigSNP"
#' @description LD pruning for a \code{bigSNP}.
#' @name Prune
#' @inheritParams bigsnpr-package
#' @param S Numeric vector of summary statistics computed
#' only with `ind.train`.
#' @param pS Numeric vector of p-values associated with S.
#' **`pS` needs to be computed through a decreasing function of `S`**.
#' For example, if `S` follows the standard normal distribution,
#' you should use `abs(S)` instead and compute
#' `pS = 2 * pnorm(abs(S), lower.tail = FALSE)`.
#' @param size Radius of the window's size for the LD evaluations.
#' @param thr.pvalue Threshold on \eqn{-log_{10}(p-value)} to assess
#' which SNPs are kept. Here, it has the purpose to accelerate computations.
#' Default is 1.
#' @param thr.corr Threshold on the correlation between two SNPs.
#' SNPs which are too correlated with another SNP which is more correlated
#' with the disease are pruned.
#' @example examples/example.pruning.R
#' @export
Prune <- function(x,
                  ind.train = seq(nrow(X)),
                  fun.stats,
                  thr.pvalue = 1,
                  size = 2000,
                  thr.corr = 0.2,
                  ncores = 1) {
  check_x(x, check.y = TRUE)

  # get descriptors
  X <- x$genotypes
  y <- x$fam$pheno
  X.desc <- describe(X)

  tmp <- fun.stats(x, ind.train)
  S <- tmp$S
  lpS <- -log10(tmp$pS)

  PruneChr <- function(lims) {
    X.chr <- sub.big.matrix(X.desc,
                            firstCol = lims[1],
                            lastCol = lims[2],
                            backingpath = x$backingpath)

    ind.chr <- seq2(lims)
    S.chr <- S[ind.chr]
    ind.col.chr <- which(lpS[ind.chr] > thr.pvalue)
    ind.keep <- list()
    l <- Inf
    while (l > 0) {
      ind <- ind.col.chr[which.max(S.chr[ind.col.chr])]
      ind.keep[length(ind.keep) + 1] <- ind.chr[ind]

      ind.col.chr.tmp <- intersect(ind.col.chr, ind + -size:size)

      res <- R_squared_chr(pBigMat = X.chr@address,
                           rowInd = ind.train,
                           colInd = ind.col.chr.tmp,
                           colMat0 = X.chr[, ind])

      ind.del <- ind.col.chr.tmp[res > thr.corr]
      ind.col.chr <- setdiff(ind.col.chr, ind.del)
      l <- length(ind.col.chr)
      #print(l)
    }

    sort(unlist(ind.keep))
  }

  range.chr <- LimsChr(x)

  obj <- foreach::foreach(i = 1:nrow(range.chr),
                          .combine = 'c',
                          .noexport = 'x')
  expr_fun <- function(i) PruneChr(range.chr[i, ])
  foreach2(obj, expr_fun, ncores)
}
