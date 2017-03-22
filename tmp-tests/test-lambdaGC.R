require(bigsnpr)

test <- snp_attachExtdata()
G <- test$genotypes
gwas <- big_univLogReg(G, sample(0:1, nrow(G), TRUE))
snp_manhattan(gwas, map = test$map)
snp_qq(gwas)

xtr <- attr(gwas, "transfo")(gwas$score)
f.opt <- function(x) (attr(gwas, "predict")(x) - 0.5)^2
lamGC <- median(xtr) / optimize(f.opt, interval = range(xtr))$minimum

lamGC
snp_qq(gwas, main = "Q-Q plot")
legend("topleft", x.intersp=0, legend = eval(substitute(expression(lambda[GC] == l),
                                           list(l = lamGC))))
plot(1:10)

snp_qq(snp_gc(gwas), main = "Q-Q plot")
