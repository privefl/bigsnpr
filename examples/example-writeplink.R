N <- 17
M <- 911

fake <- snp_fake(N, M)
fake$genotypes[] <- sample(c(0:2, NA), size = N * M, replace = TRUE)

# write the object as a bed/bim/fam object
tmpfile <- tempfile()
bed <- snp_writeBed(fake, paste0(tmpfile, ".bed"))
# read this new file for the first time
fake2 <- snp_attach(snp_readBed(bed, backingfile = basename(tmpfile),
                                backingpath = dirname(tmpfile)))

# same content
print(all.equal(fake$genotypes[,], fake2$genotypes[,]))
print(all.equal(fake$fam, fake2$fam))
print(all.equal(fake$map, fake2$map))

# two different files
print(fake$backingfile != fake2$backingfile)
