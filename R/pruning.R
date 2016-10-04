#' @title LD pruning for a "bigSNP"
#' @description LD pruning for a \code{bigSNP}.
#' @name Prune
#' @inheritParams bigsnpr-package
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
                  ind.train = NULL,
                  size = 2000,
                  thr.pvalue = 1,
                  thr.corr = 0.2,
                  ncores = 1) {
  if (class(x) != "bigSNP") stop("x must be a bigSNP")

  # get descriptors
  X <- x$genotypes
  y <- x$fam$pheno
  if (is.null(y)) stop("Please use \"GetPhenos\" to get phenotypes right.")
  X.desc <- describe(X)

  if (is.null(ind.train)) ind.train <- 1:nrow(X)

  R2 <- bigstatsr::RsqClass(X, y, ind.train)
  lpS <- -log10(stats::pchisq(length(ind.train) * R2, 1, lower.tail = F))

  PruneChr <- function(lims) {
    X.chr <- sub.big.matrix(X.desc,
                            firstCol = lims[1],
                            lastCol = lims[2],
                            backingpath = x$backingpath)

    ind.chr <- seq2(lims)
    R2.chr <- R2[ind.chr]
    ind.col <- which(lpS[ind.chr] > thr.pvalue)
    ind.col.chr <- which(ind.chr %in% ind.col)
    ind.keep <- list()
    l <- Inf
    while (l > 0) {
      ind <- ind.col.chr[which.max(R2.chr[ind.col.chr])]
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
