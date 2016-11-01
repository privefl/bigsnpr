require(bigsnpr)

celiac <- AttachBigSNP("celiac")

print(system.time(
  test <- BigToBed(celiac, "test.bed")
))
# laptop: 21 sec
# 2 sec for reading/writing
