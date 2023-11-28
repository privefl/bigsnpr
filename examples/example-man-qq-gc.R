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
