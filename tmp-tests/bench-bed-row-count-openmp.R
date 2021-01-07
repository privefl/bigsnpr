library(bigsnpr)
obj.bed <- bed("../Dubois2010_data/FinnuncorrNLITUK3hap550_QC.bed")
system.time(
  count1 <- bed_counts(obj.bed, byrow = TRUE)
)

system.time(
  count2 <- bed_counts(obj.bed, byrow = TRUE, ncores = 4)
)

Rcpp::sourceCpp('src/bed-fun.cpp')
system.time(
  count3 <- bed_row_counts_cpp(obj.bed, rows_along(obj.bed), cols_along(obj.bed),
                               ncores = 4)
)

system.time(
  count4 <- bed_col_counts_cpp(obj.bed, rows_along(obj.bed), cols_along(obj.bed),
                               ncores = 4)
)
