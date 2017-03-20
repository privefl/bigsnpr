require(bigsnpr)

test <- snp_attachExtdata()
G <- test$genotypes
gwas <- big_univLogReg(G, sample(0:1, nrow(G), TRUE))
manhattan2(gwas, map = test$map)
qq2(gwas)

xtr <- attr(gwas, "transfo")(gwas$z.score)
f.opt <- function(x) (attr(gwas, "predict")(x) - 0.5)^2
lamGC <- median(xtr) / optimize(f.opt, interval = range(xtr))$minimum

lamGC
legend("topleft", legend = eval(substitute(expression(lambda[GC] == l),
                                           list(l = lamGC))))
