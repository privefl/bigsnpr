require(bigsnpr)

bedfile <- "../thesis-celiac/Dubois2010_data/FinnuncorrNLITUK3hap550.bed"

print(system.time(
  test <- snp_readBed(bedfile, backingfile = "testRead")
)) # 19 sec -> 7.5 sec

test <- snp_attach("backingfiles/testRead.rds")

print(system.time(
  test2 <- snp_writeBed(test, "backingfiles/tmp/test5.bed")
)) # 4.7 sec

print(system.time(
  test3 <- snp_attach(snp_readBed(test2, backingfile = "testRead2"))
))

all.equal(test, test3) # filename and savedIn

all.equal(attach.BM(test$genotypes)[, 1:1e4],
          attach.BM(test3$genotypes)[, 1:1e4])
