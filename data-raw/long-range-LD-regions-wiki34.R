url <- "http://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_%28LD%29"
tables <- XML::readHTMLTable(url, header = TRUE, which = 1,
                             colClasses = c(rep("integer", 3), "character"),
                             stringsAsFactors = FALSE)
names(tables) <- c("Chr", "Start", "Stop", "ID")
LD.wiki34 <- tables

devtools::use_data(LD.wiki34)
