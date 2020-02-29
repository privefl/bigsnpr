library(bigsnpr)

celiac <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")
G <- celiac$genotypes

ind_chr <- which(celiac$map$chromosome == 2)

system.time(
  corr <- snp_cor(G, ind.col = ind_chr,
                  infos.pos = celiac$map$physical.pos[ind_chr])
) # 32 sec

system.time(
  corr2 <- snp_cor(G, ind.col = ind_chr, ncores = 4,
                   infos.pos = celiac$map$physical.pos[ind_chr])
) # 18 sec

all.equal(corr, corr2)
