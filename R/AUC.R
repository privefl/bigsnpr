################################################################################

#' AUC
#'
#' Compute the Area Under the ROC Curve (AUC) of a predictor
#' and possibly its 95\% confidence interval.
#'
#' @details Other packages provide ways to compute the AUC.
#' I chose to compute the AUC through its statistical definition as a
#' probability: \deqn{\P{score(x_{case}) > score(x_{control})}}.
#' Note that I consider equality between scores as a 50%-probability of
#' one being greater than the other.
#'
#' @return The AUC, a probability, and possibly its 2.5\% and 97.5\% quantiles
#' (95\% CI).
#'
#' @param pred Vector of predictions.
#' @param target Vector of true labels (must have exactly two levels,
#' no missing values).
#' @param nboot Number of bootstrap samples to evaluate the 95\% CI.
#' @param nsim Number of pairs (case/control) that are compared
#' in order to approximate the probability.
#' @param seed See \code{\link{set.seed}}. Use it for reproducibility.
#' @param digits See \code{\link{round}}.
#' @seealso [pROC::auc] \cr [AUC::auc]
#' @references Tom Fawcett. 2006. An introduction to ROC analysis.
#' Pattern Recogn. Lett. 27, 8 (June 2006), 861-874.
#' \url{http://dx.doi.org/10.1016/j.patrec.2005.10.010}.
#' @examples
#' snp_aucSample(c(0, 0), 0:1) # Equality of scores
#' snp_aucSample(c(0.2, 0.1, 1), c(-1, -1, 1)) # Perfect AUC
#' x <- rnorm(100)
#' z <- rnorm(length(x), x, abs(x))
#' y <- sign(z)
#' print(snp_aucSample(x, y))
#' print(snp_aucSampleConf(x, y, nboot = 1e4, nsim = 1e3))
#' @name auc
NULL

################################################################################

#' @name snp_aucSample
#' @rdname auc
#' @export
snp_aucSample <- function(pred, target, nsim = 1e7, seed = NULL, digits = 3) {
  y <- transform_levels(target)

  pred.case    <- pred[y == 1]
  pred.control <- pred[y == 0]

  set.seed(seed)

  pred.case.sample    <- sample(pred.case,    nsim, replace = TRUE)
  pred.control.sample <- sample(pred.control, nsim, replace = TRUE)

  sup <- mean(pred.case.sample > pred.control.sample)
  equ <- mean(pred.case.sample == pred.control.sample) / 2

  round(sup + equ, digits)
}

################################################################################

#' @name snp_aucSampleConf
#' @rdname auc
#' @export
snp_aucSampleConf <- function(pred, target, nboot = 1e4, nsim = 1e4,
                          seed = NULL, digits = 3) {
  y <- transform_levels(target)

  pred.case    <- pred[y == 1]
  pred.control <- pred[y == 0]

  set.seed(seed)

  sampleRes <- function() {
    ind.case <- sample(length(pred.case), replace = TRUE)
    ind.control <- sample(length(pred.control), replace = TRUE)

    pred.case.sample    <- sample(pred.case[ind.case],
                                  nsim, replace = TRUE)
    pred.control.sample <- sample(pred.control[ind.control],
                                  nsim, replace = TRUE)

    sup <- mean(pred.case.sample > pred.control.sample)
    equ <- mean(pred.case.sample == pred.control.sample) / 2

    sup + equ
  }

  res <- replicate(nboot, sampleRes())
  q <- stats::quantile(res, c(0.025, 0.975))

  round(c("Mean" = mean(res), q, "Sd" = sd(res)), digits)
}

################################################################################
