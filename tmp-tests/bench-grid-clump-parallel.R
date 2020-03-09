library(bigsnpr)

celiac <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")
G <- celiac$genotypes$copy(code = c(0, 1, 2, 0, rep(NA, 252)))
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos

gwas <- big_univLinReg(G, celiac$fam$affection)
lpval <- -predict(gwas)

system.time(
  keep <- snp_grid_clumping(G, CHR, POS, lpS = lpval,
                            exclude = which(!CHR %in% 19:22),
                            grid.base.size = 100,
                            ncores = 4)
)
# Before: 34 sec with 1 core and 20 sec with 4 cores -> parallelize chromosomes
# OpenMP: 33 sec with 1 core and 13 sec with 4 cores

# 68 sec for all 4 base size

# 140 sec for all chromosomes
