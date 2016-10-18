#' Title
#'
#' @param x
#' @param fun.model
#' @param K
#' @param progress
#' @param ...
#'
#' @return
#' @export
#' @import foreach
#'
#' @examples
crossval.ext <- function(x, fun.model, K) {
  n <- nrow(x$genotypes)
  ind <- sample(rep_len(1:K, n))

  foreach(i = 1:K) %do% {
    ind.test <- which(ind == i)
    ind.train <- setdiff(1:n, ind.test)

    fun.model(x, ind.train, ind.test, fun.model = fun.model)
  }
}

crossval.int <- function(x, fun.model, ind.train, K) {
  n <- length(ind.train)
  ind <- sample(rep_len(1:K, n))

  res.all <- foreach(i = 1:K, .combine = 'cbind') %do% {
    ind.train.test <- ind.train[ind == i]
    ind.train.train <- setdiff(ind.train, ind.train.test)

    fun.model(x, ind.train.train, ind.train.test, fun.model = NULL)
  }

  which.max(rowSums(res.all))
}
