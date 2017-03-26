require(bigsnpr)
require(pcadapt)

test <- snp_attachExtdata()
svd <- big_SVD(test$genotypes, snp_scaleBinom(), k = 10)
tmp <- snp_pcadapt(test$genotypes, svd$u[, 1:3])
plot(tmp) + scale_y_log10(limits = c(1, 12))

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

x <- tmp
m <- length(lpval <- -log10(predict(x)))
coeff <- 1

lamGC <- signif(getLambdaGC(tmp), 4)
expr <- substitute(expression(lambda[GC] == l), list(l = lamGC))
require(ggplot2)
MY_THEME(qplot(x = -log10(ppoints(m)), y = sort(lpval, decreasing = TRUE)),
          title = "Q-Q plot", coeff = coeff) +
  xlab(XLAB) + ylab(YLAB) +
  geom_abline(slope = 1, intercept = 0, color = "red", lwd = 1) +
  labs(subtitle = eval(expr))

qplot(x = lpval, geom = "qq")

plot(1:10, main = YLAB)


print(p <- snp_qq(tmp))
p$layers[[1]] <- NULL
p2 <- p + geom_point(aes(text = asPlotlyText(test$map)))
plotly::ggplotly(p2, tooltip = "text")
