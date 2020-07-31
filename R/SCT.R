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
#' @param groups List of vectors of indices to define your own categories.
#'   This could be used e.g. to derive C+T scores using two different GWAS
#'   summary statistics, or to include other information such as functional
#'   annotations. Default just makes one group with all variants.
#' @param exclude Vector of SNP indices to exclude anyway.
#'
#' @import Matrix
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
  groups = list(cols_along(G)),
  exclude = NULL,
  ncores = 1
) {

  check_args()
  assert_lengths(cols_along(G), infos.chr, infos.pos, infos.imp, lpS)
  assert_class(groups, "list")

  THR_IMP        <- sort(unique(grid.thr.imp))
  THR_CLMP       <- sort(unique(grid.thr.r2))
  BASE_SIZE_CLMP <- sort(unique(grid.base.size))

  grid <- expand.grid(
    size    = BASE_SIZE_CLMP,
    thr.r2  = THR_CLMP,
    grp.num = seq_along(groups),
    thr.imp = THR_IMP
  )
  grid$size <- as.integer(grid$size / grid$thr.r2)

  ind.noexcl <- setdiff(seq_along(infos.chr), exclude)

  all_keep <- lapply(split(ind.noexcl, infos.chr[ind.noexcl]), function(ind.chr) {

    ind.keep <- list()

    # Init chromosome
    info.chr  <- infos.imp[ind.chr]
    S.chr     <- lpS[ind.chr]
    pos.chr   <- infos.pos[ind.chr]
    stats     <- snp_colstats(G, ind.row, ind.chr, ncores)
    sumX.chr  <- stats$sumX
    denoX.chr <- stats$denoX
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

      for (group in groups) {

        ind2 <- which(ind.chr %in% group)

        # TODO: get rid of empty groups? Would be harder to compute C+T and SCT
        if (length(ind2) == 0) {

          for (thr.clmp in THR_CLMP) {
            for (base.size.clmp in BASE_SIZE_CLMP) {
              ind.keep[[length(ind.keep) + 1L]] <- integer()
            }
          }

        } else {

          ind.sp.grp   <- ind2 - 1L
          ind.chr.grp  <- ind.chr[ind2]
          ord.chr.grp  <- order(S.chr[ind2], decreasing = TRUE)
          rank.chr.grp <- match(seq_along(ord.chr.grp), ord.chr.grp)
          pos.chr.grp  <- pos.chr[ind2]
          sumX.chr.grp <- sumX.chr[ind2]
          deno.chr.grp <- denoX.chr[ind2]

          keep <- FBM(1, length(ind.chr.grp), type = "integer")

          # Loop over parameters of clumping
          for (thr.clmp in THR_CLMP) {
            for (base.size.clmp in BASE_SIZE_CLMP) {

              keep[] <- -1  # reinit to "do not know yet"

              spcor.chr <- clumping_chr_cached(
                G,
                keep,
                spcor.chr,
                spInd   = ind.sp.grp,
                rowInd  = ind.row,
                colInd  = ind.chr.grp,
                ordInd  = ord.chr.grp,
                rankInd = rank.chr.grp,
                pos     = pos.chr.grp,
                sumX    = sumX.chr.grp,
                denoX   = deno.chr.grp,
                size    = 1000 * base.size.clmp / thr.clmp, # in bp
                thr     = thr.clmp,
                ncores  = ncores
              )

              stopifnot(all(keep[] %in% 0:1))

              ind.keep[[length(ind.keep) + 1L]] <- ind.chr.grp[keep[] == 1]
            }
          }
        }
      }
    }

    ind.keep
  })

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
  exp(
    seq(from = log(from), to = log(to), length.out = length.out)
  )
}

################################################################################

#' Grid of PRS
#'
#' Polygenic Risk Scores for a grid of clumping and thresholding parameters.
#'
#' @param all_keep Output of `snp_grid_clumping()` (indices passing clumping).
#' @param betas Numeric vector of weights (effect sizes from GWAS) associated
#'   with each variant (column of `G`). If alleles are reversed, make sure to
#'   multiply corresponding effects by `-1`.
#' @param lpS Numeric vector of `-log10(p-value)` associated with `betas`.
#' @param n_thr_lpS Length for default `grid.lpS.thr`. Default is `50`.
#' @param grid.lpS.thr Sequence of thresholds to apply on `lpS`.
#'   Default is a grid (of length `n_thr_lpS`) evenly spaced on a logarithmic
#'   scale, i.e. on a log-log scale for p-values.
#' @param backingfile Prefix for backingfiles where to store scores of C+T.
#'   As we typically use a large grid, this can result in a large matrix so that
#'   we store it on disk. Default uses a temporary file.
#' @param type Type of backingfile values. Either `"float"` (the default) or
#'   `"double"`. Using `"float"` requires half disk space.
#'
#' @return `snp_grid_PRS()`: An `FBM` (matrix on disk) that stores the C+T scores
#'   for all parameters of the grid (and for each chromosome separately).
#'   It also stores as attributes the input parameters `all_keep`, `betas`,
#'   `lpS` and `grid.lpS.thr` that are also needed in `snp_grid_stacking()`.
#'
#' @rdname SCT
#' @export
snp_grid_PRS <- function(
  G, all_keep, betas, lpS,
  n_thr_lpS = 50,
  grid.lpS.thr = 0.9999 * seq_log(max(0.1, min(lpS, na.rm = TRUE)),
                                  max(lpS, na.rm = TRUE), n_thr_lpS),
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

  bigparallelr::register_parallel(ncores)

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

  res <- structure(
    multi_PRS,
    lpS = lpS,
    grid.lpS.thr = grid.lpS.thr,
    betas = betas,
    all_keep = all_keep
  )

  saveRDS(res, multi_PRS$rds)
  res
}

################################################################################

#' Stacking over PRS grid
#'
#' Stacking over many Polygenic Risk Scores, corresponding to a grid of many
#' different parameters for clumping and thresholding.
#'
#' @param multi_PRS Output of `snp_grid_PRS()`.
#' @param y.train Vector of phenotypes. If there are two levels (binary 0/1),
#'   it uses [big_spLogReg()] for stacking, otherwise [big_spLinReg()].
#' @param alphas Vector of values for grid-search. See [big_spLogReg()].
#'   Default for this function is `c(1, 0.01, 0.0001)`.
#' @param ... Other parameters to be passed to [big_spLogReg()]. For example,
#'   using `covar.train`, you can add covariates in the model with all C+T scores.
#'   You can also use `pf.covar` if you do not want to penalize these covariates.
#'
#' @rdname SCT
#' @export
snp_grid_stacking <- function(multi_PRS, y.train,
                              alphas = c(1, 0.01, 0.0001),
                              ncores = 1, ...) {

  lpS       <- attr(multi_PRS, "lpS")
  lpS_thr   <- attr(multi_PRS, "grid.lpS.thr")
  beta_gwas <- attr(multi_PRS, "betas")
  all_keep  <- attr(multi_PRS, "all_keep")

  suppressWarnings(
    mod <- `if`(length(unique(y.train)) == 2, big_spLogReg, big_spLinReg)(
      multi_PRS, y.train, alphas = alphas, ncores = ncores, ...
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
