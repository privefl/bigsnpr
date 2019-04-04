################################################################################

#' Stacked C+T (SCT)
#'
#' @inheritParams bigsnpr-package
#'
#' @name SCT
NULL

################################################################################

#' Grid of clumping
#'
#' @param grid.thr.r2 Grid of thresholds over the squared correlation between
#'   two SNPs for clumping. Default is `c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95)`.
#' @param grid.base.size Grid for base window sizes. Sizes are then computed as
#'   `base.size / thr.r2` (in kb). Default is `c(50, 100, 200, 500)`.
#' @param infos.imp Vector of imputation scores. Default is all `1` if you do
#'   not provide it.
#' @param grid.thr.imp Grid of thresholds over `infos.imp` (default is `1`), but
#'   you should change it (e.g. `c(0.3, 0.6, 0.9, 0.95)`) if providing `infos.imp`.
#' @param exclude Vector of SNP indices to exclude anyway.
#'
#' @rdname SCT
#' @export
snp_grid_clumping <- function(
  G, infos.chr, infos.pos, lpS,
  ind.row = rows_along(G),
  grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
  grid.base.size = c(50, 100, 200, 500),
  infos.imp = rep(1, ncol(G)),
  grid.thr.imp = 1,
  exclude = NULL,
  ncores = 1
) {

  check_args()
  assert_lengths(cols_along(G), infos.chr, infos.pos, infos.imp)

  THR_IMP        <- sort(unique(grid.thr.imp))
  THR_CLMP       <- sort(unique(grid.thr.r2))
  BASE_SIZE_CLMP <- sort(unique(grid.base.size))
  ALL_CHR        <- sort(unique(infos.chr))

  if (ncores == 1) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  all_keep <- foreach(chr = ALL_CHR, .packages = "Matrix") %dopar% {

    ind.keep <- list(); i <- 1

    # Init chromosome
    ind.chr   <- setdiff(which(infos.chr == chr), exclude)
    info.chr  <- infos.imp[ind.chr]
    S.chr     <- lpS[ind.chr]
    pos.chr   <- infos.pos[ind.chr]
    stats     <- big_colstats(G, ind.row = ind.row, ind.col = ind.chr)
    sumX.chr  <- stats$sum
    denoX.chr <- (length(ind.row) - 1) * stats$var
    spcor.chr <- sparseMatrix(i = integer(), j = integer(), x = double(),
                              dims = rep(length(ind.chr), 2))
    assert_sorted(pos.chr)

    for (thr.imp in THR_IMP) {

      # Init threshold of imputation
      ind <- which(info.chr >= thr.imp)

      ind.chr   <- ind.chr[ind]
      spcor.chr <- spcor.chr[ind, ind, drop = FALSE]
      info.chr  <- info.chr[ind]
      pos.chr   <- pos.chr[ind]
      S.chr     <- S.chr[ind]
      sumX.chr  <- sumX.chr[ind]
      denoX.chr <- denoX.chr[ind]

      # Loop over parameters of clumping
      for (thr.clmp in THR_CLMP) {
        for (base.size.clmp in BASE_SIZE_CLMP) {

          res <- clumping_chr_cached(
            G, spcor.chr,
            rowInd = ind.row,
            colInd = ind.chr,
            ordInd = order(S.chr, decreasing = TRUE),
            pos    = pos.chr,
            sumX   = sumX.chr,
            denoX  = denoX.chr,
            size   = 1000 * base.size.clmp / thr.clmp, # in bp
            thr    = thr.clmp
          )

          ind.keep[[i]] <- ind.chr[res[[1]]]
          i <- i + 1
          spcor.chr <- res[[2]]
        }
      }
    }

    ind.keep
  }

  grid <- expand.grid(
    size    = BASE_SIZE_CLMP,
    thr.r2  = THR_CLMP,
    thr.imp = THR_IMP
  )
  grid$size <- as.integer(grid$size / grid$thr.r2)

  structure(all_keep, grid = grid)
}

################################################################################

#' Sequence, evenly spaced on a logarithmic scale
#'
#' @inheritParams base::seq
#'
#' @examples
#' seq_log(1, 1000, 4)
#' seq_log(1, 100, 5)
#'
#' @export
seq_log <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}

#' Grid of PRS
#'
#' Polygenic Risk Scores for a grid of clumping and thresholding parameters.
#'
#' @param all_keep Output of `snp_grid_clumping()` (indices passing clumping).
#' @param betas Numeric vector of weights (effect sizes from GWAS) associated
#'   with each variant (column of `G`). If alleles are reversed, make sure to
#'   multiply corresponding effects by `-1`.
#' @param lpS Numeric vector of `-log10(p-value)` associated with `betas`.
#' @param grid.lpS.thr Sequence of thresholds to apply on `lpS`.
#'   You may want to make a grid evenly spaced on a logarithmic scale,
#'   i.e. on a log-log scale for p-values.
#' @param backingfile Prefix for backingfiles where to store scores of C+T.
#'   As we typically use a large grid, this can result in a large matrix so that
#'   we store it on disk. Default uses a temporary file.
#' @param type Type of backingfile values. Either `"float"` (the default) or
#'   `"double"`. Using `"float"` requires half disk space.
#'
#' @rdname SCT
#' @export
snp_grid_PRS <- function(
  G, all_keep, betas, lpS,
  grid.lpS.thr = seq_log(0.1, 0.999 * max(lpS), 50),
  ind.row = rows_along(G),
  backingfile = tempfile(),
  type = c("float", "double"),
  ncores = 1
) {

  check_args()
  assert_lengths(cols_along(G), betas, lpS)

  all_keep2 <- unlist(all_keep, recursive = FALSE)
  n_thr <- length(grid.lpS.thr)
  multi_PRS <- FBM(nrow = length(ind.row),
                   ncol = length(all_keep2) * n_thr,
                   backingfile = backingfile,
                   type = match.arg(type))

  if (ncores == 1) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  foreach(ic = seq_along(all_keep2)) %dopar% {
    ind.keep <- all_keep2[[ic]]
    prs <- snp_PRS(
      G, ind.test = ind.row, ind.keep = ind.keep,
      betas.keep = betas[ind.keep],
      lpS.keep = lpS[ind.keep], thr.list = grid.lpS.thr
    )
    without_downcast_warning(
      multi_PRS[, (ic - 1L) * n_thr + seq_len(n_thr)] <- prs)
    NULL
  }

  structure(
    multi_PRS$save(),
    lpS = lpS,
    grid.lpS.thr = grid.lpS.thr,
    betas = betas,
    all_keep2 = all_keep2
  )
}

################################################################################

#' Stacking over PRS grid
#'
#' Stacking over many Polygenic Risk Scores, corresponding to a grid of many
#' different parameters for clumping and thresholding.
#'
#' @param multi_PRS Output of `snp_grid_PRS()`. It stores the C+T scores for
#'   all parameters of the grid, and the `rds` file for accessing C+T scores
#'   stored on disk (that also have the dimension of chromosomes).
#'   It also stores as attributes the input parameters `all_keep`, `betas`,
#'   `lpS` and `grid.lpS.thr` that are also needed in this function.
#' @param y.train Vector of phenotypes. If there are two levels (binary 0/1),
#'   it uses [big_spLogReg()] for stacking, otherwise [big_spLinReg()].
#' @param covar.train Matrix of covariates.
#' @param pf.covar A multiplicative factor for the penalty applied to each
#'   covariate. Default does not penalize covariates (factors equal to `0`).
#' @param alphas Vector of values for grid-search. See [big_spLogReg()].
#' @param ... Other parameters to be passed to [big_spLogReg()].
#'
#' @rdname SCT
#' @export
snp_grid_stacking <- function(multi_PRS, y.train,
                              covar.train = matrix(0, length(y.train), 0L),
                              pf.covar = rep(0, ncol(covar.train)),
                              alphas = 10^(-(0:4)),
                              ncores = 1, ...) {

  lpS       <- attr(multi_PRS, "lpS")
  lpS_thr   <- attr(multi_PRS, "grid.lpS.thr")
  beta_gwas <- attr(multi_PRS, "betas")
  all_keep2 <- attr(multi_PRS, "all_keep2")

  suppressWarnings(
    mod <- `if`(length(unique(y.train)) == 2, big_spLogReg, big_spLinReg)(
      multi_PRS, y.train, alphas = alphas, ncores = ncores,
      covar.train = covar.train, pf.covar = pf.covar, ...
    )
  )

  best_mod <- summary(mod, best.only = TRUE)
  beta_best_mod <- best_mod$beta[[1]]
  beta_stacking <- rep(0, ncol(multi_PRS))
  ind_col <- attr(mod, "ind.col")
  beta_stacking[ind_col] <- head(beta_best_mod, length(ind_col))

  ind_last_thr <- 1L + sapply(lpS, function(lp) sum(lp > lpS_thr))
  coef <- rep(0, length(beta_gwas))
  n_thr_pval <- length(lpS_thr)
  ind <- seq_len(n_thr_pval)
  for (ind.keep in all_keep2) {
    b <- beta_stacking[ind]
    b2 <- c(0, cumsum(b))
    coef[ind.keep] <- coef[ind.keep] + b2[ind_last_thr[ind.keep]]
    ind <- ind + n_thr_pval
  }

  list(
    intercept  = best_mod$intercept,
    beta.G     = coef * beta_gwas,
    beta.covar = tail(beta_best_mod, -length(ind_col)),
    mod = mod
  )
}

################################################################################
