#' @name R2all
#' @title Marginally compute R2 for each column.
#' @export
R2all <- function(X, y, ind.train = seq(nrow(X)), weighted = FALSE) {
  if (weighted) {
    prop.case <- mean(y[ind.train] == 1)
    ratio <- 1 / prop.case - 1
    weights <- ifelse(y[ind.train] > 0, ratio, 1)
    R2 <- R_squared(X@address, y, ind.train, weights)
  } else {
    R2 <- R_squared(X@address, y, ind.train, rep(1, length(ind.train)))
  }

  return(R2)
}



#' @name betasall
#' @title Marginally compute betas of linear regression $y ~ X$ for each column.
#' @export
betasall <- function(X, y, ind.train = seq(nrow(X)), weighted = FALSE) {
  if (weighted) {
    prop.case <- mean(y[ind.train] == 1)
    ratio <- 1 / prop.case - 1
    weights <- ifelse(y[ind.train] > 0, ratio, 1)
    R2 <- betasRegLin(X@address, y, ind.train, weights)
  } else {
    R2 <- betasRegLin(X@address, y, ind.train, rep(1, length(ind.train)))
  }

  return(R2)
}

