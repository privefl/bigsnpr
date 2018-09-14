library(bigsnpr)

bed <- "../Dubois2010_data/FinnuncorrNLITUK3hap550.bed"
system.time(
  test <- snp_readBed(bed, tempfile())
)
# 54 / 78
G <- snp_attach(test)$genotypes
system.time(big_counts(G))
system.time(big_counts(G))
system.time(big_counts(G))

n <- fpeek::peek_count_lines(sub("\\.bed$", ".fam", bed))
p <- fpeek::peek_count_lines(sub("\\.bed$", ".bim", bed))
c(n, p)

system.time(FBM.code256(n, p, code = bigsnpr:::CODE_012, init = 2L))
# 35 sec

system.time({
  X <- FBM.code256(n, p, code = bigsnpr:::CODE_012)
  system.time(
    read_plink(X, bed, n, p, decode = bigsnpr:::getCode())
  )
})
# 41 / 54

# system.time({
#   X <- FBM.code256(n, p, code = bigsnpr:::CODE_012)
#   big_parallelize(X, p.FUN = function(X, ind, read, bed, n, p) {
#     read(X, bed, n, p, decode = bigsnpr:::getCode(), ind_col = ind)
#   }, p.combine = unlist, ind = seq_len(p), ncores = 2,
#   read = bigsnpr:::read_plink, bed = bed, n = n, p = p)
# })
# 49 / 65
# parallelism is useless, as always when writing to disk

