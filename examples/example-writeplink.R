N <- 17
M <- 911

fake <- snp_fake(N, M)
X <- attach.BM(fake$genotypes)
X[] <- sample(as.raw(0:3), length(X), TRUE)
rm(X)

# write the object as a bed/bim/fam object
tmp <- tempfile(fileext = ".bed")
bed <- snp_writeBed(fake, tmp)
# read this new file for the first time
rds <- snp_readBed(bed, backingfile = "test_write")
# attach object in R session
fake2 <- snp_attach(rds)

# same content
all.equal(attach.BM(fake$genotypes)[,],
          attach.BM(fake2$genotypes)[,])
all.equal(fake$fam, fake2$fam)
all.equal(fake$map, fake2$map)

# two different files
fake$savedIn
fake2$savedIn



# cleaning
bk <- sub("\\.rds$", ".bk", rds)
file.remove(rds, bk)
