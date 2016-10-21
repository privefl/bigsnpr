################################################################################

#' Create a function with specified parameters.
#'
#' @param fun The function you want to use.
#' @param ... Replace this by a list of parameters of
#' `fun` you want to set.
#'
#' @return A new function, derived from `fun`.
#' @export
#' @seealso [crossval.ext] [crossval.int]
#' @examples
#' round3 <- model.fun(round, digits = 3)
#' round3(4.5588877)
model.fun <- function(fun, ...) {
  args.fun <- list(...)
  function(...) do.call(fun, args = c(list(...), args.fun))
}

################################################################################

#' K-fold cross-validation.
#'
#' K-fold cross-validation of `fun.model` on `x`.
#'
#' @inheritParams bigsnpr-package
#' @param fun.model Function that is used to train and
#' predict a particular model. Its first 3 parameters
#' should be `x`, `ind.train` and `ind.test` and it should
#' have a parameter `fun.model` with `NULL` as default value.
#' @param K Number of turns in the cross-validation.
#'
#' @return A list of each turn's score.
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

################################################################################

#' K-fold cross-validation.
#'
#' Internal K-fold cross-validation of `fun.model` on `x`
#' in order to tune hyper-parameters.
#'
#' @inheritParams crossval.ext
#' @param ind.train Indices used as training set by the model.
#'
#' @return The position of the parameter which has given
#' the best scores (in average) of the model.
#' @export
#' @import foreach
#'
#' @examples
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

################################################################################
