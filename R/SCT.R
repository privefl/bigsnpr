################################################################################

#' Stacked C+T (SCT)
#'
#' @name SCT
NULL

################################################################################

#' Grid of clumping
#'
#' @rdname SCT
#' @export
snp_grid_clumping <- function(G, infos.chr, infos.pos, lpS,
                              ind.row = rows_along(G),
                              grid.thr.r2 = c(0.01, 0.05, 0.2, 0.8, 0.95),
                              grid.base.size = c(30, 100, 300),
                              infos.imp = NULL,
                              grid.thr.imp = NULL,
                              exclude = NULL,
                              ncores = 1) {

  check_args()

  if (is.null(infos.imp)) infos.imp <- rep(1, ncol(G))
  if (is.null(grid.thr.imp)) grid.thr.imp <- 1
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

    for (thr.imp in THR_IMP) {

      # Init threshold of imputation
      ind <- which(info.chr >= thr.imp)

      ind.chr   <- ind.chr[ind]
      spcor.chr <- spcor.chr[ind, ind]
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
    thr.imp = THR_IMP,
    chr     = ALL_CHR
  )
  grid$size <- as.integer(grid$size / grid$thr.r2)

  structure(all_keep, grid = grid)
}

################################################################################

#' Grid of PRS
#'
#' Polygenic Risk Scores for a grid of clumping and thresholding parameters.
#'
#' @rdname SCT
#' @export
snp_grid_PRS <- function(
  G, all_keep, betas, lpS,
  grid.lpS.thr = exp(seq(log(0.1), log(0.999 * max(lpS)), length.out = 50)),
  ind.row = rows_along(G),
  same_ref = rep(TRUE, length(betas)),
  backingfile = tempfile(),
  ncores = 1
) {

  check_args()
  assert_lengths(cols_along(G), betas, lpS, same_ref)

  n_thr_pval <- length(grid.lpS.thr)
  grid_size <- length(all_keep[[1]]) * n_thr_pval
  scores_all_chr <- matrix(0, length(ind.row), grid_size)
  scores_by_chr <- FBM(length(ind.row), length(all_keep) * grid_size,
                       backingfile = backingfile)$save()

  if (ncores == 1) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  offset <- 0L
  for (ind_keep in all_keep) {
    all_prs_chr <- foreach(ind.keep = ind_keep, .combine = "cbind") %dopar% {
      snp_PRS(
        G, ind.test = ind.row, ind.keep = ind.keep,
        betas.keep = betas[ind.keep], same.keep = same_ref[ind.keep],
        lpS.keep = lpS[ind.keep], thr.list = grid.lpS.thr
      )
    }
    scores_all_chr <- scores_all_chr + all_prs_chr
    scores_by_chr[, offset + seq_len(grid_size)] <- all_prs_chr
    offset <- offset + grid_size
  }

  structure(
    scores_all_chr,
    rds = scores_by_chr$rds,
    lpS = lpS,
    grid.lpS.thr = grid.lpS.thr,
    betas = betas,
    all_keep = all_keep
  )
}

################################################################################

#' Stacking over PRS grid
#'
#' @rdname SCT
#' @export
snp_grid_stacking <- function(multi_PRS, y.train, all_keep,
                              covar.train = matrix(0, length(y.train), 0L),
                              pf.covar = rep(0, ncol(covar.train)),
                              alphas = 10^(-(0:4)),
                              ncores = 1) {

  rds       <- attr(multi_PRS, "rds")
  lpS       <- attr(multi_PRS, "lpS")
  lpS_thr   <- attr(multi_PRS, "grid.lpS.thr")
  beta_gwas <- attr(multi_PRS, "betas")
  all_keep  <- attr(multi_PRS, "all_keep")

  # scores_all_chr <- multi_PRS
  scores_by_chr <- big_attach(rds)

  suppressWarnings(
    mod <- `if`(length(unique(y.train)) == 2, big_spLogReg, big_spLinReg)(
      scores_by_chr, y.train, alphas = alphas, ncores = ncores,
      covar.train = covar.train, pf.covar = pf.covar
    )
  )

  best_mod <- summary(mod, best.only = TRUE)
  beta_best_mod <- best_mod$beta[[1]]
  beta_stacking <- rep(0, ncol(scores_by_chr))
  ind_col <- attr(mod, "ind.col")
  beta_stacking[ind_col] <- head(beta_best_mod, length(ind_col))

  ind_last_thr <- 1L + sapply(lpS, function(lp) sum(lp > lpS_thr))
  coef <- rep(0, length(beta_stacking))
  n_thr_pval <- length(lpS_thr)
  ind <- seq_len(n_thr_pval)
  for (ind.keep in unlist(all_keep, recursive = FALSE)) {
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
