################################################################################

#' Manhattan plot
#'
#' Creates a manhattan plot.
#'
#' If you don't have information of chromosome and position, you should simply
#' use `plot` instead.
#'
#' @param gwas A `mhtest` object with the p-values associated with each SNP.
#' Typically, the output of [bigstatsr::big_univLinReg], [bigstatsr::big_univLogReg]
#' or [snp_pcadapt].
#' @inheritParams bigsnpr-package
#' @param colors Colors used for each chromosome (they are recycled).
#' Default is an alternation of black and gray.
#' @param dist.sep.chrs "Physical" distance that separates two chromosomes.
#' Default is 10 Mbp.
#' @param ind.highlight Indices of SNPs you want to highlight (of interest).
#' Default doesn't highlight any SNPs.
#' @param col.highlight Color used for highlighting SNPs. Default uses red.
#' @param npoints Number of points to keep (ranked by p-value) in order to get
#' a lighter object (and plot). Default doesn't cut anything.
#' If used, the resulting object will have an attribute called `subset`
#' giving the indices of the kept points.
#' @param labels Labels of the x axis. Default uses the number of the
#' chromosome there are in `infos.chr`(`sort(unique(infos.chr))`). This may be
#' useful to restrict the number of labels so that they are not overlapping.
#' @param coeff Relative size of text. Default is `1`.
#'
#' @return A `ggplot2` object. You can plot it using the `print` method.
#' You can modify it as you wish by adding layers. You might want to read
#' [this chapter](https://r4ds.had.co.nz/data-visualisation.html)
#' to get more familiar with the package **ggplot2**.
#' @export
#' @import ggplot2
#'
#' @example examples/example-man-qq-gc.R
snp_manhattan <- function(gwas, infos.chr, infos.pos,
                          colors = c("black", "grey60"),
                          dist.sep.chrs = 1e7,
                          ind.highlight = integer(0),
                          col.highlight = "red",
                          labels = NULL,
                          npoints = NULL,
                          coeff = 1) {

  check_args()
  assert_lengths(infos.chr, infos.pos, rows_along(gwas))

  # get ordering
  ord <- order(infos.chr, infos.pos)
  infos.chr <- infos.chr[ord]
  infos.pos <- infos.pos[ord]

  # get all chromosomes
  all.chr <- sort(unique(infos.chr))
  if (is.null(labels)) labels <- all.chr

  # get plotting positions of each SNP and chromosome
  offset <- 0
  all.pos <- lapply(all.chr, function(chr) {
    ind.chr <- which(infos.chr == chr)
    pos <- infos.pos[ind.chr] + offset + dist.sep.chrs
    offset <<- tail(pos, 1)
    pos
  })
  label.pos <- sapply(all.pos, function(pos) mean(range(pos)))
  all.pos <- unlist(all.pos)

  # get colors
  colors <- rep_len(colors, length(all.chr))
  all.colors <- colors[match(infos.chr, all.chr)]
  all.colors[ind.highlight] <- col.highlight

  # get plot
  lpval <- stats::predict(gwas)[ord]
  cond <- is.null(npoints)
  ind <- if (cond) seq_along(lpval) else head(order(lpval), npoints)
  ymin <- -lpval[tail(ind, 1)]
  subtitle <- paste0("(values \u2265 ", signif(ymin), ")")  # \u2265: >=
  p <- ggplot(data.frame(pos = all.pos, lp = -lpval)[ind, ], aes(pos, lp)) +
    geom_point(color = all.colors[ind]) +
    scale_x_continuous(breaks = label.pos, labels = labels,
                       limits = range(all.pos)) +
    labs(title = "Manhattan Plot", x = "Chromosome",
         y = expression(-log[10](italic("p-value"))),
         subtitle = `if`(cond, NULL, subtitle)) +
    theme_bigstatsr(size.rel = coeff)

  p$plot_env <- emptyenv()

  `if`(cond, p, structure(p, subset = ind))
}

################################################################################

getLambdaGC <- function(gwas, tol = 1e-8) {

  check_args()

  xtr <- attr(gwas, "transfo")(stats::na.omit(gwas)[["score"]])
  PREDICT <-  attr(gwas, "predict")
  MEDIAN <- log10(0.5)
  f.opt <- function(x) (PREDICT(x) - MEDIAN)

  stats::median(xtr) / stats::uniroot(f.opt, interval = range(xtr),
                                      check.conv = TRUE, tol = tol)$root
}

################################################################################

#' Q-Q plot
#'
#' Creates a quantile-quantile plot from p-values from a GWAS study.
#'
#' @param lambdaGC Add the Genomic Control coefficient as subtitle to the plot?
#'
#' @inherit snp_manhattan return params
#' @export
#' @import ggplot2
#'
#' @example examples/example-man-qq-gc.R
snp_qq <- function(gwas, lambdaGC = TRUE, coeff = 1) {

  check_args()

  p <- graphics::plot(gwas, type = "Q-Q", coeff = coeff)

  if (lambdaGC) {
    lamGC <- signif(getLambdaGC(gwas))
    expr <- substitute(expression(lambda[GC] == l), list(l = lamGC))
    p + labs(subtitle = eval(expr))
  } else {
    p
  }
}

################################################################################

#' Genomic Control
#'
#' @export
#'
#' @inherit snp_manhattan return params
#'
#' @references Devlin, B., & Roeder, K. (1999).
#' Genomic control for association studies.
#' Biometrics, 55(4), 997-1004.
#'
#' @example examples/example-man-qq-gc.R
snp_gc <- function(gwas) {

  check_args()

  force(lamGC <- getLambdaGC(gwas))

  # https://stackoverflow.com/a/42938212/6103040
  gcf <- function(f, lamGC) {
    transfo <- f
    function(x) transfo(x) / lamGC
  }
  attr(gwas, "transfo") <- gcf(attr(gwas, "transfo"), lamGC)

  gwas
}

################################################################################
