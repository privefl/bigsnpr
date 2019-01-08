library(bigsnpr)
prefix <- "tmp-data/indep"
plink <- download_plink()
bed2 <- snp_plinkQC(plink, prefix.in = prefix,
                    maf = 0.05, geno = 0.05, mind = 0.05, autosome.only = TRUE)
bed3 <- snp_plinkIBDQC(plink, bed2, ncores = nb_cores(), pi.hat = 0.12,
                       pruning.args = c(500, 0.1))
rds <- snp_readBed(bed3)
snp <- snp_attach(rds)
G <- snp$genotypes
CHR <- snp$map$chromosome
counts <- big_counts(G)
counts[, 1:10]
ind <- which(counts[4, ] > 0 & counts[4, ] <= 16)
big_apply(G, a.FUN = function(X, ind) {
  block <- X[, ind]
  matrixStats::colMedians(G[, 1:10], na.rm = TRUE)
})
infos <- snp_fastImpute(G, CHR, seed = 1, ncores = nb_cores())
# In the meantime, in another session, you can check progress
infos <- big_attach("tmp-data/indep-infos-impute.rds")
mean(!is.na(infos[1, ]))
