#' @name PCA.bigSNP
#' @title Principal Components of a "bigSNP".
#' @description Get k or all Principal Components (PCs) of a \code{bigSNP}
#' @inheritParams bigstatsr::BigXYt
#' @param x A \code{bigSNP}.
#' @param k Number of PCs to compute. Default is all.
#' @param thr.eigval Threshold to remove "unsignificant" PCs.
#' Default is \code{1e-3}.
#' @param scaling Scaling to be used, either one of
#' \itemize{
#' \item "none",
#' \item "center": center columns,
#' \item "scale": center columns and divide by sd,
#' \item "binom": use \eqn{(x-2p)/sqrt(2p(1-p))} (Patterson 2006),
#' where p is the MAF,
#' \item "y-aware": center columns and  multiply by beta of
#' univariate linear regression.
#' See \href{https://goo.gl/8G8WMa}{this blog post}.
#' }
#' @references Patterson N, Price AL, Reich D (2006) Population Structure and
#' Eigenanalysis. PLoS Genet 2(12): e190. doi: 10.1371/journal.pgen.0020190
#' @export
#' @return A \code{matrix} of PCs.
#' @details See \code{\link[bigstatsr]{BigXYt}}.
#'
#' Note that for the Eigen decomposition, only \code{R} is
#' used because is faster (see \href{http://goo.gl/UYJcCw}{stackoverflow}).
#' If you want a large number of eigenvectors/values, you should
#' really considerer using Microsoft R Open for speed.
#' @example examples/example.PCA.bigSNP.R
#' @seealso \code{\link{bigSNP}} \code{\link{prcomp}}
PCA.bigSNP <- function(x,
                       block.size,
                       k = NULL,
                       ind.train = seq(nrow(X)),
                       scaling = "none",
                       thr.eigval = 1e-3,
                       use.Eigen = TRUE,
                       progress = TRUE) {
  check_x(x)
  X <- x$genotypes

  if (scaling == "none") {
    vec.center <- rep(0, length(ind.train))
    vec.scale <- rep(1, length(ind.train))
  } else if (scaling == "center") {
    vec.center <- bigstatsr::colmeans(X, ind.train)
    vec.scale <- rep(1, length(ind.train))
  } else if (scaling == "scale") {
    tmp <- bigstatsr::colmeans_sds(X, ind.train)
    vec.center <-  tmp$mean
    vec.scale <- tmp$sd
  } else if (scaling == "binom") {
    p <- bigstatsr::colmeans(X, ind.train) / 2
    vec.center <- 2*p
    vec.scale <- sqrt(2*p*(1 - p))
  } else if (scaling == "y-aware") {
    vec.center <- bigstatsr::colmeans(X, ind.train)
    y <- x$pheno
    if (all(sort(unique(y)) == c(-1, 1))) {
      vec.scale <- 1 / bigstatsr::CoeffsClass(X, y, ind.train)
    } else {
      vec.scale <- 1 / bigstatsr::CoeffsReg(X, y, ind.train)
    }
  } else {
    stop("This is not implemented yet.")
  }

  rotated <- bigstatsr::BigPCA(X = X,
                               block.size = block.size,
                               k = k,
                               ind.train = ind.train,
                               vec.center = vec.center,
                               vec.scale = vec.scale,
                               thr.eigval = thr.eigval,
                               use.Eigen = use.Eigen,
                               progress = progress)

  return(rotated)
}

