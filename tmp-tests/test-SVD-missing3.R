library(bigsnpr)
popres <- snp_attach("../paper-packages/backingfiles/POPRESQC.rds")

G <- popres$genotypes
CHR <- popres$map$chromosome
POS <- popres$map$physical.pos

ind.keep <- snp_clumping(G, infos.chr = CHR, infos.pos = POS, ncores = 2,
                         exclude = snp_indLRLDR(CHR, POS))
obj.svd <- big_randomSVD(G, snp_scaleBinom(), ind.col = ind.keep,
                         ncores = 2)

# Compute minor allele frequencies
maf <- snp_MAF(G)

# Generate number of missing values by SNP with a Beta-Binomial distribution
n <- nrow(G)
m <- ncol(G)
nbNA <- VGAM::rbetabinom.ab(m, size = n, shape1 = 0.6, shape2 = 5)
sum(nbNA) / (n * m)

# Generate indices of missing values
indNA <- cbind(
  unlist(lapply(nbNA, function(nb) {
    `if`(nb > 0, sample(n, size = nb), NULL)
  })),
  rep(cols_along(G), nbNA)
)

# Fill a copy of the matrix with NAs (coded as 03)
GNA <- big_copy(G)
GNA[indNA] <- as.raw(3)

# Write a new bedfile
popresNA <- popres
popresNA$genotypes <- GNA
bedfile <- snp_writeBed(popresNA, ind.col = ind.keep,
                        bedfile = tempfile(fileext = ".bed"))

obj.bed <- bed(bedfile)
obj.svd2 <- bed_randomSVD(obj.bed)
obj.svd3 <- bed_randomSVD(obj.bed, ncores = 2)
library(pcadapt)
obj.svd4 <- pcadapt(read.pcadapt(bedfile, type = "bed"), K = 10, pca.only = TRUE)

# GNA$backingfile
# infos <- snp_fastImpute(GNA, infos.chr = CHR, p.train = 0.9, ncores = 2)
# obj.svd5 <- big_randomSVD(GNA$copy(CODE_IMPUTE_PRED), snp_scaleBinom(),
#                           ind.col = ind.keep, ncores = 2)

# In another session:
# library(bigstatsr)
# back <- "/tmp/RtmpycNRlX/file51e231ccaeeb"
# system(paste0("ls ", back, "*"))
# infos <- big_attach(paste0(back, "-infos-impute.rds"))
# infos[, 1:5]
# # df <- data.frame(pNA = infos[1, ], pError = infos[2, ])
#
# repeat {
#   p <- mean(!is.na(infos[1, ]))
#   cat(round(100 * p, 1), "- ")
#   if (p == 1) break else Sys.sleep(60)
# }

rbind(obj.svd2$d - obj.svd$d,
      obj.svd3$d - obj.svd$d,
      obj.svd4$d - obj.svd$d,
      obj.svd5$d - obj.svd$d)
# [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
# [1,]  11.49680  16.548036  20.329226  21.256263  22.061032  23.556689
# [2,] -65.98938 -39.962300 -30.566604 -27.949647 -26.258241 -24.887522
# [3,]  10.78069  15.964330  19.824155  20.887771  21.677397  23.159696
# [4,] -10.46815  -6.079897  -3.797146  -2.957417  -2.310035  -1.655515
# [,7]       [,8]       [,9]      [,10]
# [1,]  23.462814  23.217800  23.379923  23.321418
# [2,] -24.656254 -24.814958 -24.640958 -24.568134
# [3,]  23.061166  22.838750  22.997012  22.911263
# [4,]  -1.858613  -2.069427  -1.775503  -1.719268

c(mean(sqrt(colSums(cor(obj.svd2$u, obj.svd$u)^2))),
  mean(sqrt(colSums(cor(obj.svd3$u, obj.svd$u)^2))),
  mean(sqrt(colSums(cor(obj.svd4$u, obj.svd$u)^2))),
  mean(sqrt(colSums(cor(obj.svd5$u, obj.svd$u)^2))))
# 0.8271456 0.8471444 0.8246797 0.8937024
# 0.8641700 0.8748350 0.8632781 0.9110531
