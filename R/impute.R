################################################################################

#' Imputation
#'
#' Imputation via local XGBoost models.
#'
#' This implementation aims at imputing with good accuracy while being fast.
#' For instance, `nrounds` and `max_depth` are voluntarily chosen low.
#' One could use for example `nrounds = 50` and `max_depth = 5` to improve
#' accuracy a bit, at the expense of computation time.
#' \cr
#' Moreover, `breaks` and `sizes` are used to dynamically choose the number
#' of predictors to use more complex models when there are more missing values.
#' To use only complex models, one could use `breaks = 0` and `sizes = 100`,
#' increasing computation time.
#'
#' @inheritParams bigsnpr-package
#' @param perc.train Percentage of individuals used for the training
#' (the rest is used to assess the error of imputation by SNP).
#' @param nrounds The number of iterations (trees) in
#' [xgboost][xgboost::xgboost]. Default is `20`.
#' @param max_depth maximum depth of a tree. Default is `3`.
#' @param breaks Break points of the number of missing values by SNPs
#' to which we change the number of predictors to be used. The more missing
#' values there are, the more complex is the model used.
#' @param sizes Radius of how many predictors enter the model.
#' @param verbose Print progress? Default is `FALSE`.
#'
#' @return The new [bigSNP] object with a slot `imputation` which is
#' a `data.frame` with 2 columns:
#' - the number of missing values by SNP,
#' - the estimated error of imputation by SNP.
#' @export
#' @import xgboost
#'
#' @examples
snp_impute <- function(x, perc.train = 0.7, nrounds = 20, max_depth = 3,
                       breaks = c(0, trunc(n / c(1000, 200, 100))),
                       sizes = c(10, 20, 30, 50),
                       ncores = 1, verbose = FALSE) {
  check_x(x)
  X <- x$genotypes
  n <- nrow(X)
  n.train <- round(n * perc.train)
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
    error <- numeric(m.part)
    # many different `ind.train` and `ind.val`
    n.rep <- 2 * max(sizes)
    ind.train.rep <- replicate(n.rep, sort(sample(n, n.train)))
    ind.val.rep <- apply(ind.train.rep, 2, function(ind) setdiff(seq(n), ind))

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
      # if (!(i %% 10)) print(i)
      X.label <- X.part[, i] * 1
      nbNA[i] <- l <- length(indNA <- which(is.na(X.label)))
      if (l > 0) {
        j <- i %% n.rep + 1
        ind <- setdiff(ind.train.rep[, j], indNA)
        ind.val <- ind.val.rep[, j]

        s <- sizes[max(which(l > breaks))]
        X.data <- X.part[, interval(i, s), drop = FALSE] * 1

        bst <- xgboost(data = X.data[ind, , drop = FALSE],
                       label = X.label[ind],
                       nrounds = nrounds, params = params,
                       nthread = 1, verbose = 0)

        pred <- predict(bst, X.data[indNA, , drop = FALSE])
        X2.part[indNA, i] <- round2(pred)
        pred2 <- predict(bst, X.data[ind.val, , drop = FALSE])
        error[i] <- mean(round2(pred2) != X.label[ind.val], na.rm = TRUE)
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
