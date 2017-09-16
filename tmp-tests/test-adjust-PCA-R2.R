# Doesn't work
library(bigsnpr)

rds <- snp_readBed("../POPRES_data/POPRES_allchr_QC_norel.bed",
                   "backingfiles/POPRESQC")
popres <- snp_attach(rds)
G <- popres$genotypes
CHR <- popres$map$chromosome
POS <- popres$map$physical.pos

cor <- snp_split(popres$map$chromosome,
                 function(G, ind.chr) {
                   cor <- snp_cor(G, ind.col = ind.chr)
                 }, ncores = 2, G = G)


ind.chr2 <- which(CHR == 2)
cor <- snp_cor(G, ind.col = ind.chr2)
library(Matrix)
r2.sqrt <- colSums(cor^2)
plot(r2.sqrt, pch = 19, cex = .6,
     col = 1 + seq_along(tmp) %in% snp_indLRLDR(CHR[ind.chr2], POS[ind.chr2]))

test <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, ind.col = ind.chr2)
plot(test)
plot(test, type = "scores")

test_scale <- function(X, ind.row, ind.col) {
  tmp <- snp_scaleBinom()(X, ind.row, ind.col)
  tmp$scale <- tmp$scale * r2.sqrt
  tmp
}

test2 <- big_SVD(G, test_scale, ind.col = ind.chr2)
plot(test2)
plot(test2, type = "scores")
plot(test2, type = "scores", scores = 3:4)
