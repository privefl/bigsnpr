library(bigsnpr)

celiac <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")
G <- celiac$genotypes

ind_chr <- which(celiac$map$chromosome == 22)

GNA <- big_copy(G, ind.col = ind_chr[1:200])
set.seed(1); ind <- sort(sample(length(GNA), length(GNA) / 100))
GNA[ind] <- 3
infos.imp <- bigsnpr:::FBM_infos(GNA)
infos.imp$rds

# system.time(
#   infos <- snp_fastImpute(GNA, rep(1, ncol(GNA)))
# ) # 183 sec

system.time(
  bigsnpr:::imputeChr(GNA, infos.imp, cols_along(GNA),
                      alpha = 1e-4,
                      size = 200,
                      p.train = 0.8,
                      n.cor = nrow(GNA),
                      seed = NA,
                      ncores = 4)
)
## N = 5000 / M = 2000
# 1 -> 178 sec / 4 -> 114 sec
## N = All / M = 1000
# 1 -> 253 sec / 4 -> 138 sec / 6 -> 138 sec

# 1350 / 26 -> 1350 / 13


GNA <- big_copy(G)
set.seed(1); ind <- sort(sample(length(GNA), length(GNA) / 100)); GNA[ind] <- 3
system.time(
  snp_fastImputeSimple(GNA, ncores = 4)
)
