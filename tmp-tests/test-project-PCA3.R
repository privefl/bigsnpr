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
attr(test, "shrinkage")
obj.svd.ref <- attr(test, "obj.svd.ref")
PC.ref <- predict(obj.svd.ref)
points(PC.ref, col = "red", pch = 20)
plot(test[, 5:6])
points(PC.ref[, 5:6], col = "red", pch = 20)
plot(test[, 5:6 + 2])
points(PC.ref[, 5:6 + 2], col = "red", pch = 20)

bedfile.new2 <- "../POPRES_data/POPRES_allchr.bed"
test2 <- bed_projectPCA(
  bed.new = bed(bedfile.new2),
  bed.ref = bed(bedfile.ref),
  join_by_pos = TRUE,
  build.new = "hg18", liftOver = "tmp-data/liftOver",
  ncores = nb_cores()
)

plot(test2)
str(attributes(test2))
attr(test2, "shrinkage")
obj.svd.ref2 <- attr(test2, "obj.svd.ref")
PC.ref2 <- predict(obj.svd.ref2)
points(PC.ref2, col = "red", pch = 20)
plot(test2[, 5:6])
points(PC.ref2[, 5:6], col = "red", pch = 20)
plot(test2[, 5:6 + 2])
points(PC.ref2[, 5:6 + 2], col = "red", pch = 20)

library(ggplot2)
dist <- bigsnpr:::robLOF(test2, c(4, 10, 30), ncores = nb_cores())
qplot(test2[, 1], test2[, 2], color = dist) +
  scale_colour_viridis_c() +
  theme_bigstatsr()
# LOF should be preferred if want to keep many pops
# RobMaha should be preferred if want to keep a single population
hist(dist2 <- sqrt(dist), "FD")
abline(v = print(q <- bigutilsr::tukey_mc_up(dist2)), col = "red")
qplot(test2[, 1], test2[, 2], color = (dist2 > q)) +
  scale_colour_viridis_d() +
  theme_bigstatsr()
