wd <- setwd("~/Bureau/bigsnpr/tmp-data")

library(bigsnpr)

bedfile.ref <- download_1000G(".")
bed.ref <- bed(bedfile.ref)
ref.map <- setNames(bed.ref$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
ref.map2 <- snp_modifyBuild(ref.map, "./liftOver", from = "hg19", to = "hg18")

bed.new <- bed("../../POPRES_data/POPRES_allchr.bed")
new.map <- setNames(bed.new$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
bigsnpr::snp_match(cbind(ref.map2, beta = 1), new.map, join_by_pos = FALSE)
nrow(.Last.value)  # 0
bigsnpr::snp_match(cbind(ref.map2, beta = 1), new.map)
nrow(.Last.value)  # 304,951

new.map2 <- snp_modifyBuild(new.map, "./liftOver", from = "hg18", to = "hg19")
bigsnpr::snp_match(cbind(ref.map, beta = 1), new.map2)
nrow(.Last.value)  # 304,953

setwd(wd)
