LD.wiki34 <- XML::readHTMLTable(
  "http://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_%28LD%29",
  header = TRUE,
  which = 1,
  colClasses = c(rep("integer", 3), "character"),
  stringsAsFactors = FALSE
)
names(LD.wiki34) <- c("Chr", "Start", "Stop", "ID")

# devtools::use_data(LD.wiki34)
