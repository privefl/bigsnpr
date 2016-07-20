################################################################################

#' @title AUC.
#' @description Compute the Area Under the ROC Curve (AUC) of a predictor
#' and possibly its 95\% confidence interval.
#' @details Other packages provide ways to compute the AUC.
#' Yet, I found out that when there are equalities, those can
#' overestimate the true AUC. So, I choose to compute the AUC
#' through an (accurate) estimation of its statistical definition as a
#' probability: \deqn{\P{score(x_{case}) > score(x_{control})}}
#' @return The AUC, a probability, and possibly its 2.5\% and 97.5\% quantiles
#' (95\% CI).
#' @param pred Vector of predictions.
#' @param y Vector of true values (for the same individuals).
#' @param nboot Number of bootstrap samples to evaluate 95\% CI.
#' @param nsim Number of pairs (case/control) that are compared
#' in order to approximate the probability.
#' @param seed See \code{\link{set.seed}}. Use it for reproducibility.
#' @param digits See \code{\link{round}}.
#' @seealso pROC::\code{\link[pROC]{auc}} \cr
#' AUC::\code{\link[AUC]{auc}} \cr
#' Hmisc::\code{\link[Hmisc]{somers2}}
#' @references Tom Fawcett. 2006. An introduction to ROC analysis.
#' Pattern Recogn. Lett. 27, 8 (June 2006), 861-874.
#' DOI=\url{http://dx.doi.org/10.1016/j.patrec.2005.10.010}.
#' @name auc
NULL

################################################################################

#' @name AucSample
#' @rdname auc
#' @export
AucSample <- function(pred, y, nsim = 1e7, seed = NULL, digits = 3) {
  pred.case    <- pred[y == 1]
  pred.control <- pred[y == -1]

  set.seed(seed)

  pred.case.sample    <- sample(pred.case,    nsim, replace = TRUE)
  pred.control.sample <- sample(pred.control, nsim, replace = TRUE)

  return(round(mean(pred.case.sample > pred.control.sample), digits))
}

################################################################################

#' @name AucSampleConf
#' @rdname auc
#' @export
AucSampleConf <- function(pred, y, nboot = 1e5, nsim = 1e4,
                         seed = NULL, digits = 3) {
  pred.case    <- pred[y == 1]
  pred.control <- pred[y == -1]

  set.seed(seed)

  sampleRes <- function() {
    ind.case <- sample(length(pred.case), replace = TRUE)
    ind.control <- sample(length(pred.control), replace = TRUE)
    pred.case.sample    <- sample(pred.case[ind.case],
                                  nsim, replace = TRUE)
    pred.control.sample <- sample(pred.control[ind.control],
                                  nsim, replace = TRUE)
    mean(pred.case.sample > pred.control.sample)
  }

  res <- replicate(nboot, sampleRes())
  q <- stats::quantile(res, c(0.025, 0.975))
  return(round(c("Mean" = mean(res), q), digits))
}

################################################################################
