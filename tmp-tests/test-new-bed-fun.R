bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
# bedfile <- "../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel.bed"
bimfile <- sub("\\.bed$", ".bim", bedfile)
famfile <- sub("\\.bed$", ".fam", bedfile)
n_total <- bigreadr::nlines(famfile)
m_total <- bigreadr::nlines(bimfile)
ind.row <- seq_len(n_total)
chr <- bigreadr::fread2(bimfile, select = 1)[[1]]
ind.col <- which(chr == 1)
lookup_byte <- bigsnpr:::getCode()
Rcpp::sourceCpp('src/bed-fun.cpp')
stats <- bedcolvars(bedfile, n_total, m_total, ind.row, ind.col, lookup_byte)
center <- stats$sum / stats$nona
scale <- sqrt(stats$var * (stats$nona - 1))
lookup_scale <- rbind(t(sapply(0:2, function(g) (g - center) / scale)), 0)
ord <- sample(length(ind.col))
pos <- bigreadr::fread2(bimfile, select = 4)[[1]]
keep <- bed_clumping_chr(bedfile, n_total, m_total, ind.row, ind.col,
                         lookup_byte, lookup_scale,
                         ord, pos[ind.col], 500e3, 0.2)
mean(keep)

