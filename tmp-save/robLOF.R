################################################################################

utils::globalVariables("nPC")

#' @examples
#' X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
#' pca <- prcomp(X, scale. = TRUE, rank. = 10)
#' U <- pca$x
#' robLOF(U)
#' library(ggplot2)
#' qplot(U[, 1], U[, 2], color = robLOF(U)) + coord_equal() +
#'   scale_color_viridis_c()
robLOF <- function(U, seq_kNN,
                   combine = max,
                   robMaha = FALSE,
                   log = TRUE,
                   ncores = 1) {

  if (ncores == 1) {
    registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  all_llof <- foreach(nPC = rev(tail(cols_along(U), -1)), .combine = "cbind") %dopar% {
    bigutilsr::LOF(U[, 1:nPC], seq_k = seq_kNN, robMaha = robMaha, log = log,
                   combine = combine)
  }

  scale <- apply(all_llof, 2, bigutilsr::tukey_mc_up)
  all_llof_scaled <- sweep(all_llof, 2, scale, '/')
  apply(all_llof_scaled, 1, combine)
}

################################################################################
