library(bigsnpr)

bedfile <- "../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.bed"

gc()
system.time({
  rds <- snp_readBed(bedfile, tempfile())
  print(rds)
}, gcFirst = FALSE) # 12-13 sec

gc()
system.time({
  rds <- snp_readBed2(bedfile, tempfile())
  print(rds)
}, gcFirst = FALSE) # 18 sec

gc()
system.time({
  rds <- snp_readBed2(bedfile, tempfile(), ncores = 4)
  print(rds)
}, gcFirst = FALSE) # 8 sec

# ind.col = rep(1:200e3, each = 5) -> 27 sec / 31 sec with collapse(2)
