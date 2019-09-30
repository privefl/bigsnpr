library(bigsnpr)
obj.bed <- bed("../POPRES_data/POPRES_allchr.bed"); nPC <- 10
obj.svd <- bed_autoSVD2(obj.bed, k = nPC, ncores = nb_cores())
plot(obj.svd)
attr(obj.svd, "lrldr")
plot(obj.svd, type = "loadings", loadings = 3:8, coef = 0.5)
plot(obj.svd, type = "scores", scores = 3:10, coef = 0.5)
