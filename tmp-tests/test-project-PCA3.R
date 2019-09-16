bedfile.new <- "../Dubois2010_data/FinnuncorrNLITUK3hap550.bed"
bedfile.ref <- download_1000G("tmp-data")

test <- bed_projectPCA(
  bed.new = bed(bedfile.new),
  bed.ref = bed(bedfile.ref),
  join_by_pos = FALSE,
  ncores = nb_cores()
)

plot(test)
str(attributes(test))

bedfile.new2 <- "../POPRES_data/POPRES_allchr.bed"
test2 <- bed_projectPCA(
  bed.new = bed(bedfile.new2),
  bed.ref = bed(bedfile.ref),
  join_by_pos = TRUE,
  build.new = "hg18", liftOver = "tmp-data/liftOver",
  ncores = nb_cores()
)
str(attributes(test2))
plot(test2[, 1:2])
plot(test2[, 3:4])
plot(test2[, 5:6])
plot(test2[, 7:8])
plot(test2[, 9:10])

library(ggplot2)
dist <- bigutilsr::LOF(test2[, 1:2], 30)
qplot(test2[, 1], test2[, 2], color = dist) +
  scale_colour_viridis_c() +
  theme_bigstatsr()
# LOF should be preferred if want to keep many pops
# RobMaha should be preferred if want to keep a single population
