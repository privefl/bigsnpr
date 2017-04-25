set.seed(1)

test <- snp_attachExtdata()
G <- test$genotypes
n <- nrow(X)
m <- ncol(X)

# pruning / clumping with MAF
ind.keep <- snp_pruning(G, infos.chr = test$map$chromosome, thr.r2 = 0.1)
# keep most of them -> not much LD in this simulated dataset
length(ind.keep) / m

ind.keep2 <- snp_clumping(G, infos.chr = test$map$chromosome, thr.r2 = 0.1)
