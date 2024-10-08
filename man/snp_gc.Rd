% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/man-qq-gc.R
\name{snp_gc}
\alias{snp_gc}
\title{Genomic Control}
\usage{
snp_gc(gwas)
}
\arguments{
\item{gwas}{A \code{mhtest} object with the p-values associated with each SNP.
Typically, the output of \link[bigstatsr:big_univLinReg]{bigstatsr::big_univLinReg}, \link[bigstatsr:big_univLogReg]{bigstatsr::big_univLogReg}
or \link{snp_pcadapt}.}
}
\value{
A \code{ggplot2} object. You can plot it using the \code{print} method.
You can modify it as you wish by adding layers. You might want to read
\href{https://r4ds.had.co.nz/data-visualisation.html}{this chapter}
to get more familiar with the package \strong{ggplot2}.
}
\description{
Genomic Control
}
\examples{
set.seed(9)

test <- snp_attachExtdata()
G <- test$genotypes
y <- rnorm(nrow(G))

gwas <- big_univLinReg(G, y)

snp_qq(gwas)
gwas_gc <- snp_gc(gwas) # this modifies `attr(gwas_gc, "transfo")`
snp_qq(gwas_gc)

# The next plot should be prettier with a real dataset
snp_manhattan(gwas_gc,
              infos.chr = test$map$chromosome,
              infos.pos = test$map$physical.pos) +
  ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = 2, color = "red")

p <- snp_qq(gwas_gc) +
  ggplot2::aes(text = asPlotlyText(test$map)) +
  ggplot2::labs(subtitle = NULL, x = "Expected -log10(p)", y = "Observed -log10(p)")
\dontrun{plotly::ggplotly(p, tooltip = "text")}
}
\references{
Devlin, B., & Roeder, K. (1999).
Genomic control for association studies.
Biometrics, 55(4), 997-1004.
}
