library(bigsnpr)
fake <- snp_fake(200, 100)
fake$genotypes[] <- sample(0:3, length(fake$genotypes), replace = TRUE,
                           prob = c(1, 1, 1, 3))
(n_bad <- sum(big_counts(fake$genotypes)[4, ] > 100))
bed <- snp_writeBed(fake, tempfile(fileext = ".bed"))
library(testthat)
expect_warning(ind.keep <- bed_clumping(bed, thr.r2 = 0.01),
               sprintf("%d variants have >50%% missing values.", n_bad))

plink <- download_plink("tmp-data")
tmp <- sub("\\.bed$", "", bed)
counts <- big_counts(fake$genotypes)
af <- drop(crossprod(counts[1:3, ], 0:2)) / colSums(counts[1:3, ]) / 2
maf <- pmin(af, 1 - af)
write.table(data.frame(SNP = fake$map$marker.ID, P = 1 - maf),
            file = paste0(tmp, ".frq"), row.names = FALSE, quote = FALSE)
# Clumping
library(glue)
system(glue("{plink} --bfile {tmp} --out {tmp}",
            " --clump {tmp}.frq",
            " --clump-p1 1 --clump-p2 1 --clump-r2 0.2"))
ind <- match(read.table(glue("{tmp}.clumped"), header = TRUE)$SNP, fake$map$marker.ID)
2 * length(intersect(ind, ind.keep)) / (length(ind) + length(ind.keep))
length(intersect(ind, ind.keep)) / length(union(ind, ind.keep))
