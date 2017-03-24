################################################################################

MY_THEME <- bigstatsr:::MY_THEME

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
snp_manhattan <- function(gwas, map, ...) {

  stopifnot(inherits(gwas, "mhtest"))
  stopifnot(nrow(gwas) == nrow(map))

  qqman::manhattan(data.frame(CHR = map[[1]],
                              SNP = map[[2]],
                              BP  = map[[4]],
                              P   = predict(gwas)),
                   ...)
}

################################################################################

getLambdaGC <- function(gwas) {

  stopifnot(inherits(gwas, "mhtest"))

  xtr <- attr(gwas, "transfo")(gwas[["score"]])
  f.opt <- function(x) (attr(gwas, "predict")(x) - 0.5)^2

  lamGC <- median(xtr) / optimize(f.opt, interval = range(xtr))$minimum
}

#' Q-Q plot
#'
#' Creates a quantile-quantile plot from p-values from a GWAS study based on
#' the package **qqman**.
#'
#' @param gwas A `mhtest` object with the p-values associated with each SNP.
#' Typically, the output of [big_univLinReg] or [big_univLogReg].
#' @param lambdaGC Add the Genomic Control cofficient to the plot?
#' @param ... Other arguments passed to plot().
#'
#' @return `NULL`. Creates a Q-Q plot.
#' @export
#'
#' @examples
snp_qq <- function(gwas, lambdaGC = TRUE, ...) {

  stopifnot(inherits(gwas, "mhtest"))
  pvalues <- predict(gwas)
  stopifnot(all((pvalues >= 0) & (pvalues <= 1)))
  qqman::qq(pvalues, ...)

  if (lambdaGC) {
    lamGC <- signif(getLambdaGC(gwas), 4)
    expr <- substitute(expression(lambda[GC] == l), list(l = lamGC))
    legend("topleft", x.intersp = 0, legend = eval(expr))
  }

  invisible()
}

################################################################################

#' Genomic Control
#'
#' @param gwas
#'
#' @return
#' @export
#'
#' @examples
snp_gc <- function(gwas) {

  stopifnot(inherits(gwas, "mhtest"))
  force(lamGC <- getLambdaGC(gwas))

  # http://stackoverflow.com/a/42938212/6103040
  gcf <- function(f, lamGC) {
    transfo <- f
    function(x) transfo(x) / lamGC
  }
  attr(gwas, "transfo") <- gcf(attr(gwas, "transfo"), lamGC)

  gwas
}

################################################################################
