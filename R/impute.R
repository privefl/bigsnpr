################################################################################

#' Title
#'
#' @param x
#' @param nrounds
#' @param max_depth
#' @param breaks
#' @param sizes
#' @param Kfolds
#' @param ncores
#' @param verbose
#'
#' @return
#' @export
#' @import xgboost
#'
#' @examples
snp_imputeCV <- function(x, nrounds = 20, max_depth = 3,
                         breaks = c(0, trunc(n / c(1000, 200, 100))),
                         sizes = c(10, 20, 30, 50),
                         Kfolds = c(2, 3, 5, 8),
                         ncores = 1, verbose = FALSE) {
  check_x(x)
  X <- x$genotypes
  n <- nrow(X)
  params <- list(max_depth = max_depth)

  # get descriptors
  X.desc <- describe(X)

  newfile <- checkFile(x, "impute")
  X2 <- deepcopy(X, type = "char",
                 backingfile = paste0(newfile, ".bk"),
                 backingpath = x$backingpath,
                 descriptorfile = paste0(newfile, ".desc"))
  X2.desc <- describe(X2)

  range.chr <- LimsChr(x)

  # function that imputes one chromosome
  if (is.seq <- (ncores == 1)) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores, outfile = `if`(verbose, "", NULL))
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  res <- foreach(ic = seq_len(nrow(range.chr)), .combine = 'rbind') %dopar% {
    lims <- range.chr[ic, ]

    if (verbose)
      printf("Imputing chromosome %d with \"XGBoost\"...\n", lims[3])

    X.part <- sub.big.matrix(X.desc, firstCol = lims[1], lastCol = lims[2])
    X2.part <- sub.big.matrix(X2.desc, firstCol = lims[1], lastCol = lims[2])

    m.part <- ncol(X.part)
    nbNA <- integer(m.part)
    error <- rep(NA_real_, m.part)

    # useful functions
    interval <- function(i, size) {
      ind <- i + -size:size
      ind[ind >= 1 & ind <= m.part & ind != i]
    }
    round2 <- function(pred) (pred > 0.5) + (pred > 1.5)

    opt.save <- options(bigmemory.typecast.warning = FALSE)
    on.exit(options(opt.save), add = TRUE)

    # imputation
    for (i in 1:m.part) {
      X.label <- X.part[, i] * 1
      nbNA[i] <- l.NA <- length(indNA <- which(is.na(X.label)))
      if (l.NA > 0) {
        w <- max(which(l.NA > breaks))
        ind.col <- interval(i, sizes[w])
        K <- Kfolds[w]

        X.data.noNA <- X.part[-indNA, ind.col, drop = FALSE] * 1
        X.label.noNA <- X.label[-indNA]
        X.data.NA <- X.part[indNA, ind.col, drop = FALSE] * 1

        l.noNA <- nrow(X.data.noNA)
        indCV <- sample(rep_len(1:K, l.noNA))
        preds.NA <- matrix(NA_real_, l.NA, K)
        err <- 0

        for (k in 1:K) {
          ind.val <- which(indCV == k)

          bst <- xgboost(data = X.data.noNA[-ind.val, , drop = FALSE],
                         label = X.label.noNA[-ind.val],
                         nrounds = nrounds,
                         params = params,
                         nthread = 1,
                         verbose = 0,
                         save_period = NULL)

          # validation error
          pred.val <- predict(bst, X.data.noNA[ind.val, , drop = FALSE])
          err <- err + sum(round2(pred.val) != X.label.noNA[ind.val])
          # one vector of predicted values for imputation
          preds.NA[, k] <- predict(bst, X.data.NA)
        }

        # impute by conbining predictions of each fold
        X2.part[indNA, i] <- round2(apply(preds.NA, 1, median))
        # validation error
        error[i] <- err / l.noNA
      }
    }

    if (verbose)
      printf("Done imputing chromosome %d.\n", lims[3])

    data.frame(nbNA, error)
  }

  # create new imputed bigSNP
  snp_list <- list(genotypes = X2,
                   fam = x$fam,
                   map = x$map,
                   imputation = res,
                   backingfile = newfile,
                   backingpath = x$backingpath)
  class(snp_list) <- "bigSNP"
  # save it
  saveRDS(snp_list, file.path(x$backingpath, paste0(newfile, ".rds")))
  # return it
  snp_list
}

################################################################################
