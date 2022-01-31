library(bigsnpr)

chr6 <- snp_attach("../Dubois2010_data/celiac_chr6.rds")
G <- chr6$genotypes

# debugonce(bigsnpr:::cor0)
profvis::profvis(
  corr <- snp_cor(G, ind.row = 1:1000, ind.col = 1:5000, ncores = 6)
)
# Before: 296 MB total: 27 for corMat, 160 for diag, 109 for sparseMatrix
# After: 66 MB

corr[1:5, 1:5]
corr
