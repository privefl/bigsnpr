#' @title Operations of linear regression
#' @description Operations of \bold{univariate} linear regression on all the columns of a big.matrix.\cr
#' In the case of classification, weights are used to make the results independent
#' from the proportion of cases / controls (see parameter \emph{weights} of \code{\link[stats]{lm}}).
#' @param X A \code{\link[bigmemory]{big.matrix}}.
#' @param y Either \itemize{
#' \item a vector of \{-1, 1\} in the case of classification (suffix Class),
#' \item a vector of more than two unique values in the case of regression (suffix Reg).
#' }
#' @param ind.train A vector of the row indices that are used.
#' @return \itemize{
#' \item \code{RsqReg} and \code{RsqClass} return a vector of
#' the coefficient of determination \eqn{R^2} of each regression,
#' \item \code{BetasReg} and \code{BetasClass} return a vector of
#' the slopes \eqn{\beta_1} of each regression.
#' }
#' @name reglin
NULL

#' @name RsqReg
#' @rdname reglin
#' @export
RsqReg <- function(X, y, ind.train = seq(nrow(X))) {
  len <- length(unique(y[ind.train]))
  if (len > 2) {
    return(R_squared(X@address, y, ind.train, rep(1, length(ind.train))))
  } else {
    stop(paste0("y has not enough unique elements.\n",
                "Try using RsqClass instead."))
  }
}

#' @name RsqClass
#' @rdname reglin
#' @export
RsqClass <- function(X, y, ind.train = seq(nrow(X))) {
  if (all(sort(unique(y)) == c(-1, 1))) {
    prop.case <- mean(y[ind.train] == 1)
    ratio <- 1 / prop.case - 1
    weights <- ifelse(y[ind.train] > 0, ratio, 1)
    return(R_squared(X@address, y, ind.train, weights))
  } else {
    stop("y should be a vector of 1 (cases) and -1 (controls).")
  }
}


#' @name BetasReg
#' @rdname reglin
#' @export
BetasReg <- function(X, y, ind.train = seq(nrow(X))) {
  len <- length(unique(y[ind.train]))
  if (len > 2) {
    return(betasRegLin(X@address, y, ind.train, rep(1, length(ind.train))))
  } else {
    stop(paste0("y has not enough unique elements.\n",
                "Try using BetasClass instead."))
  }
}

#' @name BetasClass
#' @rdname reglin
#' @export
BetasClass <- function(X, y, ind.train = seq(nrow(X))) {
  if (all(sort(unique(y)) == c(-1, 1))) {
    prop.case <- mean(y[ind.train] == 1)
    ratio <- 1 / prop.case - 1
    weights <- ifelse(y[ind.train] > 0, ratio, 1)
    return(betasRegLin(X@address, y, ind.train, weights))
  } else {
    stop("y should be a vector of 1 (cases) and -1 (controls).")
  }
}

