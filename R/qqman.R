################################################################################

#' Manhattan plot
#'
#' Creates a manhattan plot based on the package **qqman**.
#'
#' @param gwas A `mhtest` object with the p-values associated with each SNP.
#' Typically, the output of [big_univLinReg] or [big_univLogReg].
#' @param map The slot `map` of a [bigSNP][bigSNP-class].
#' @inheritDotParams qqman::manhattan -x -chr -bp -p -snp
#'
#' @return `NULL`. Creates a Manhattan plot.
#' @export
#'
#' @references Devlin, B., & Roeder, K. (1999).
#' Genomic control for association studies.
#' Biometrics, 55(4), 997-1004.
#'
#' @examples
manhattan2 <- function(gwas, map, ...) {

  stopifnot(nrow(gwas) == nrow(map))

  qqman::manhattan(cbind(gwas, map), chr = NAMES.MAP[1], bp = NAMES.MAP[4],
                   p = "p.value", snp = NAMES.MAP[2], ...)
}

################################################################################

#' Q-Q plot
#'
#' Creates a quantile-quantile plot from p-values from a GWAS study based on
#' the package **qqman**.
#'
#' @param pvalues A numeric vector of p-values or a data.frame (or just a list)
#' with a named vector called "p.value".
#' @param lambdaGC Add the Genomic Control cofficient to the plot?
#' @param ... Other arguments passed to plot().
#'
#' @return `NULL`. Creates a Q-Q plot.
#' @export
#'
#' @examples
qq2 <- function(pvalues, lambdaGC = TRUE, ...) {

  if (is.data.frame(pvalues) && (ncol(pvalues) > 1))
    pvalues <- pvalues[["p.value"]]

  stopifnot(all((pvalues >= 0) & (pvalues <= 1)))


  qqman::qq(pvalues, ...)

  if (lambdaGC) {
    xtr <- attr(gwas, "transfo")(gwas$z.score)
    f.opt <- function(x) (attr(gwas, "predict")(x) - 0.5)^2
    lamGC <- median(xtr) / optimize(f.opt, interval = range(xtr))$minimum

    legend("topleft", legend = eval(substitute(expression(lambda[GC] == l),
                                               list(l = lamGC))))
  }

  NULL
}

################################################################################
