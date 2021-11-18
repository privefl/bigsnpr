library(bigsnpr)
G <- snp_attachExtdata()$genotypes
rows <- sample(nrow(G), replace = TRUE)
cols <- sample(ncol(G), replace = TRUE)
system.time(test <- snp_ld_scores(G, rows, cols, ncores = 2))
system.time(true <- snp_ld_scores(G, rows, cols, ncores = 1))
system.time(snp_cor(G, rows, cols, ncores = 1))
system.time(snp_cor(G, rows, cols, ncores = 2))

all.equal(test, true)
plot(test, true)


bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
obj.bed <- bed(bedfile)
system.time(test <- bed_ld_scores(obj.bed, rows, cols, ncores = 2))
system.time(true <- bed_ld_scores(obj.bed, rows, cols, ncores = 1))
system.time(bed_cor(obj.bed, rows, cols, ncores = 1))
system.time(bed_cor(obj.bed, rows, cols, ncores = 2))

all.equal(test, true)
plot(test, true)
