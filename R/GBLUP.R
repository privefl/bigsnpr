#' ################################################################################
#'
#' #' @name GBLUP
#' #' @title gBLUP of a "bigSNP".
#' #' @description Genetic Best Linear Unbiased Predictor (gBLUP)
#' #'  of a \code{bigSNP}
#' #' @inheritParams BigXYt
#' #' @param thr.eigval Threshold to remove eigenvalues to get
#' #' stability in the inversion of the matrix.
#' #' Default is \code{1e-3}.
#' #' @export
#' #' @return A \code{vector} of predictions.
#' #' @details See \code{\link[bigstatsr]{BigXYt}}.
#' #'
#' #' Note that for the Eigen decomposition, only \code{R} is
#' #' used because is faster (see \href{http://goo.gl/UYJcCw}{stackoverflow}).
#' #' If you want a large number of eigenvectors/values, you should
#' #' really considerer using Microsoft R Open for speed.
#' #' @example examples/example.GBLUP.R
#' #' @seealso \code{\link{bigSNP}}
#' GBLUP <- function(x, block.size, ind.train,
#'                   thr.eigval = 1e-3,
#'                   use.Eigen = TRUE) {
#'   check_x(x)
#'
#'   res <- BigXYt(x, block.size, ind.train, use.Eigen)
#'
#'   bigK  <- res[[1]]
#'   bigK2 <- res[[2]]
#'   rm(res)
#'
#'   n <- nrow(bigK)
#'   means <- bigcolsumsDouble(bigK@address) / n
#'   symCenter(bigK@address, means, mean(means))
#'   colCenter(bigK2@address, means)
#'
#'   eig <- eigen(bigK[,], symmetric = TRUE)
#'
#'   m <- ncol(x$genotypes)
#'   lastEig <- max(which(eig$values > (thr.eigval * m)))
#'
#'   y.train <- x$fam$affection[ind.train]
#'   ind <- 1:lastEig
#'   tmp <- crossprod(eig$vectors[, ind], y.train - mean(y.train))
#'   tmp2 <- tmp / eig$values[ind]
#'   tmp3 <- eig$vectors[, ind] %*% tmp2
#'   pred <- bigK2[,] %*% tmp3 + mean(y.train)
#'
#'   return(pred)
#' }
#'
#' ################################################################################
