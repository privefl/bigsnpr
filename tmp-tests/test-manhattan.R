require(bigsnpr)

popres <- snp_attach("backingfiles/popresNA.rds")
G <- popres$genotypes
chrs <- popres$map$chromosome
pos <- popres$map$physical.pos
ind.excl <- snp_indLRLDR(chrs, pos)

ind.keep <- snp_clumping(G, chrs, thr.r2 = 0.2, exclude = ind.excl, ncores = 4)
G.svd <- big_randomSVD(G, fun.scaling = snp_scaleBinom(), ind.col = ind.keep,
                       ncores = 4)
require(robust)
test <- snp_pcadapt(G, G.svd$u[, 1:5])
# plot(test)

map <- popres$map
p <- snp_manhattan(test, map$chromosome, map$physical.pos)
p.sub <- snp_manhattan(test, map$chromosome, map$physical.pos, npoints = 10e3)
print(p.sub)
require(ggplot2)
# p2 <- p + aes(text = asPlotlyText(map))
# print(p3 <- p2 + ylim(c(2, NA)))
# p3$layers

ind <- attr(p.sub, "subset")
infos <- cbind(Index = rows_along(map), map[, c(1, 2, 4)])
p4 <- plotly::ggplotly(p.sub + aes(text = asPlotlyText(infos[ind, ])),
                       tooltip = "text")

p5 <- snp_qq(snp_gc(test)) + aes(text = asPlotlyText(infos)) + xlim(c(3, NA))
plotly::ggplotly(, tooltip = "text")
