################################################################################

test <- snp_attachExtdata()
G <- test$genotypes
y01 <- test$fam$affection - 1

################################################################################

# PCA -> covariables
obj.svd <- snp_autoSVD(G, infos.chr = test$map$chromosome,
                       infos.pos = test$map$physical.position)

# GWAS
gwas <- big_univLogReg(G, y01.train = y01, covar.train = obj.svd$u)
pval <- predict(gwas, log10 = FALSE)
pval2 <- readRDS(system.file("testdata", "pval.rds", package = "bigsnpr"))
expect_equal(pval, pval2, tolerance = 1e-4)

# clumping
ind.keep <- snp_clumping(G, infos.chr = test$map$chromosome,
                         S = abs(gwas$score),
                         size = 250, # as PLINK default
                         is.size.in.bp = TRUE,
                         infos.pos = test$map$physical.pos)
ind.keep2 <- readRDS(system.file("testdata", "clumping.rds",
                                 package = "bigsnpr"))
expect_gt(mean(ind.keep %in% ind.keep2), 0.98)

plink <- download_plink()
tmp <- tempfile()
# Compute MAFs
write.table(data.frame(SNP = test$map$marker.ID, P = pval),
            file = paste0(tmp, ".gwas"), row.names = FALSE, quote = FALSE)
snp_writeBed(test, paste0(tmp, ".bed"))
# Clumping
library(glue)
system(glue("{plink} --bfile {tmp} --out {tmp}",
            " --clump {tmp}.gwas",
            " --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --allow-no-sex"))
infos <- read.table(glue("{tmp}.clumped"), header = TRUE)
ind.keep3 <- sort(match(infos$SNP, test$map$marker.ID))

setdiff(ind.keep, ind.keep2)
setdiff(ind.keep, ind.keep3)
infos[grep(test$map$marker.ID[157], infos$SP2), ]
cor(G[, 157], G[, 210])^2
system(glue("{plink} --bfile {tmp} --out {tmp}",
            " --ld {test$map$marker.ID[157]} {test$map$marker.ID[210]}"))
cor(G[, 1], G[, 3])^2
infos[grep(test$map$marker.ID[385], infos$SP2), ]
cor(G[, 385], G[, 177])^2
infos[grep(test$map$marker.ID[4178], infos$SP2), ]
cor(G[, 4178], G[, 4166])^2
# The r2 values computed by --clump are based on maximum likelihood haplotype
# frequency estimates; you can use '--r2 dprime' to dump them all.

################################################################################
