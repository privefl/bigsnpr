require(bigsnpr)
require(pcadapt)

test <- snp_attachExtdata()
files <- system.file("extdata", paste0("cluster", 1:3), package = "bigsnpr")
test <- snp_getPops(test, pop.files = files,
                    col.sample.ID = 2,
                    col.family.ID = 3)
mat <- t(attach.BM(test$genotypes)[,])

svd <- big_SVD(test$genotypes, snp_scaleBinom(), k = 10)

x <- pcadapt(mat, K = 10)
plot(x,option="scores")
plot(svd)

require(ggplot2)
scores <- predict(svd)
p <- qplot(x = scores[, 1], y = scores[, 2])
population <- test$fam$family.ID



MYTHEME <- function(p, title = NULL, coeff = 1) {
  p + theme_bw() + ggtitle(title) +
    theme(plot.title   = element_text(size = rel(2.0 * coeff), hjust = 0.5),
          legend.title = element_text(size = rel(1.5 * coeff)),
          legend.text  = element_text(size = rel(1.2 * coeff)),
          axis.title   = element_text(size = rel(1.5 * coeff)),
          axis.text    = element_text(size = rel(1.2 * coeff)))
}

population <- test$fam$family.ID
MYTHEME(qplot(x = scores[, 1], y = scores[, 2]),
        title = "Scores of PCA",
        coeff = 1.2) + xlab("PC1") + ylab("PC2") +
  theme(legend.position = c(0.83, 0.15)) + geom_point(aes(color = population))

nval <- 97

qplot(y = svd$v[, 1])


x <- pcadapt(mat, K = 3)
require(robust)
tmp <- snp_pcadapt(test$genotypes, svd$u[, 1:3])
plot(x$stat, tmp[[1]])
plot(x$pvalues, predict(snp_gc(tmp)), log = "xy")
abline(0, 1, col = "red")

snp_qq(snp_gc(tmp))
plot(x, option = "qqplot")
plot(x, option = "manhattan")
snp_manhattan(snp_gc(tmp), test$map)
