################################################################################

#' Title
#'
#' @param x
#' @param ind.train
#' @param weighted
#'
#' @return
#' @export
#'
#' @examples
R2stats <- function(x, ind.train, weighted = FALSE) {
  check_x(x, check.y = TRUE)

  X <- x$genotypes
  y <- x$fam$pheno

  R2 <- bigstatsr::RsqClass(X, y, ind.train, weighted = weighted)
  S <- length(ind.train) * R2
  pS <- pchisq(S, 1, lower.tail = FALSE)
  list(S = S, pS = pS)
}

################################################################################

#' Title
#'
#' @param x
#' @param ind.train
#' @param weighted
#' @param fun.stats
#' @param thr.list
#' @param pruning
#' @parzm K
#' @param fun.model
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
PRS1 <- function(x, ind.train, ind.test,
                 fun.stats,
                 weighted = FALSE,
                 thr.list = 0,
                 pruning = FALSE,
                 K = 10,
                 fun.model = NULL,
                 ...) {
  check_x(x, check.y = TRUE)

  printf("Weigthed = %s\n", weighted)

  X <- x$genotypes
  y <- x$fam$pheno

  # compute slopes of univariate linear regressions
  betas <-
    bigstatsr::CoeffsClass(X, y, ind.train, weighted = weighted)["Slopes", ]

  # need lpS
  tmp <- fun.stats(x, ind.train)
  lpS <- -log10(tmp$pS)

  # pruning
  if (pruning) {
    ind.keep <- Prune(x, ind.train, fun.stats, ...)
    print(length(ind.keep))
    lpS[-ind.keep] <- 0 # artificially remove them
  }

  n.thr <- length(thr.list)

  if (is.null(fun.model)) {
    res.all <- foreach(j = 1:n.thr, .combine = 'c') %do% {
      ind.col <- which(lpS > thr.list[j])
      scores <- prs1(X@address, betas, ind.train, ind.col)
      AucSample(scores, y[ind.test])
    }
  } else {
    if (n.thr > 1) { # need tuning hyper-parameter
      which.thr <- crossval.int(x, fun.model, ind.train, K)
      thr.opt <- thr.list[which.thr]
      printf("thr.opt is %s\n", thr.opt)
    } else if (n.thr == 1) {
      thr.opt <- thr.list
    } else {
      stop("'thr.list' should be at least of length 1.")
    }
    ind.col <- which(lpS > thr.opt)
    scores <- prs1(X@address, betas, ind.test, ind.col)
    AucSample(scores, y[ind.test])
  }
}

################################################################################

#' Title
#'
#' @param x
#' @param S
#' @param lpS
#' @param ind.train
#' @param ind.test
#' @param thr.list
#' @param pruning
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
PRS <- function(x, betas, S, lpS,
                ind.train = seq(nrow(X)),
                ind.test,
                thr.list = 0,
                pruning = TRUE,
                ncores = 1,
                ...) {
  check_x(x)
  X <- x$genotypes
  y <- transform_levels(x$fam$affection)

  # pruning
  if (pruning) {
    ind.keep <- snp_clumping(x, S, ind.train, ncores = ncores, ...)
    print(length(ind.keep))
    lpS[-ind.keep] <- -1 # artificially remove them
  }

  n.thr <- length(thr.list)

  res.all <- foreach(j = 1:n.thr, .combine = 'c') %do% {
    ind.col <- which(lpS > thr.list[j])
    if (length(ind.col)) {
      scores <- prs1(X@address, betas, ind.train, ind.col)
      snp_aucSample(scores, y[ind.test])
    } else {
      NA
    }
  }
}
