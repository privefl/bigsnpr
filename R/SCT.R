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
snp_grid_clumping <- function(G, infos.chr, infos.pos,
                              ind.row = rows_along(G),
                              lpS = NULL,
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
    S.chr     <- S[ind.chr]
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
    thr.imp = THR_IMP,
    thr.r2  = THR_CLMP,
    chr     = ALL_CHR
  )
  grid$size <- as.integer(grid$size / grid.thr.r2)

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
  scores_by_chr <- FBM(nrow(G), grid_size * length(all_keep),
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

  structure(scores_all_chr, rds = scores_by_chr$rds)
}

################################################################################
