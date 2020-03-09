library(bigsnpr)

bigsnp <- snp_attach("../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.rds")

gc()
system.time({
  rds <- snp_subset(bigsnp, backingfile = tempfile())
  print(rds)
}, gcFirst = FALSE) # 7-8 sec -> 16 sec

gc()
system.time({
  rds <- snp_subset(bigsnp, backingfile = tempfile(), ncores = 4)
  print(rds)
}, gcFirst = FALSE) # -> 7-8 sec
