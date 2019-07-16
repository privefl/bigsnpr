library(bigsnpr)
celiac <- snp_attach("../paper2-PRS/backingfiles/celiacQC_sub1.rds")

map.1000G <- bigreadr::fread2(
  "tmp-data/1000G_phase3_common_hapmap.bim",
  col.names = c("chr", "rsid", "osef", "pos", "a1", "a0")
)

info_snp <- snp_match(
  cbind(map.1000G[-3], beta = 1),
  setNames(celiac$map[-3], c("chr", "rsid", "pos", "a1", "a0")),
  join_by_pos = FALSE
)
info_snp

bed.1000G <- bed("tmp-data/1000G_phase3_common_hapmap.bed")
svd.1000G <- bed_autoSVD(bed.1000G, ind.col = info_snp$`_NUM_ID_.ss`, k = 20,
                         ncores = nb_cores())

plot(svd.1000G)
plot(svd.1000G, type = "scores")
plot(svd.1000G, type = "scores", scores = 1:18, cols = 3, coeff = 0.7)
plot(svd.1000G, type = "scores", scores = 9:10 + 7)
plot(svd.1000G, type = "loadings", loadings = 1:10, coeff = 0.4)

info_sample <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
info_sample
sum(duplicated(info_sample$`Individual ID`))
info_sample <- info_sample[info_sample$`Individual ID` %in% bed.1000G$fam$sample.ID, ]

fam <- bed.1000G$fam
ped <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
lapply(ped$Siblings, function(x) {
  if (x == "0") return(NULL) else trimws(scan(text = x, what = "", sep = ","))
})
fam2 <- dplyr::left_join(fam[c(2, 5)], ped[c(1:5, 7)],
                         by = c("sample.ID" = "Individual ID", "sex" = "Gender"))
pop <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")
fam2 <- dplyr::left_join(fam2, pop[1:3], by = c("Population" = "Population Code"))
str(fam2)

adj <- matrix(0, nrow(fam), nrow(fam), dimnames = list(fam$sample.ID, fam$sample.ID))
ind_rel <- rbind(
  cbind(rows_along(fam), match(fam$paternal.ID, fam$sample.ID, nomatch = 0L)),
  cbind(rows_along(fam), match(fam$maternal.ID, fam$sample.ID, nomatch = 0L))
)
adj[ind_rel] <- 0.5
ind <- which(colSums(adj != 0) > 0)
g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE,
                                 diag = FALSE, add.colnames = NULL, add.rownames = NULL)
plot(igraph::induced_subgraph(g, vids = ind))
# https://en.wikipedia.org/wiki/Coefficient_of_relationship
rel_coef <- c(
  "child" = 0.5, "Child" = 0.5,
  "father" = 0.5, "mother" = 0.5, "mother; child" = 0.5,
  "mat grandfather" = 0.25, "mat grandmother" = 0.25,
  "pat grandmother" = 0.25, "pat grandfather" = 0.25,
  "pat grandfather; father" = 0.25, "pat grandmother; mother" = 0.25,
  "unrel" = 0, "unrels" = 0,
  "husband of Child" = 0, "wife of child" = 0, "not father" = 0
)


get_rel_ind <- function(i, var, coef) {
  x <- ped[[i, var]]
  if (x != "0") {
    IDs <- trimws(scan(text = x, what = "", sep = ",", quiet = TRUE))
    data.frame(i, IDs, coef, stringsAsFactors = FALSE)
  } else {
    NULL
  }
}
ind <- do.call(
  rbind,
  lapply(match(fam$sample.ID, ped$`Individual ID`), function(i) {
    rbind(get_rel_ind(i, "Siblings",     1 / 2),
          get_rel_ind(i, "Second Order", 1 / 4),
          get_rel_ind(i, "Third Order",  1 / 8))
  })
)
rel_to_fam <- ped$`Individual ID`[ind[ind[[2]] %in% fam$sample.ID, 1]]
is_rel <- (fam$sample.ID %in% rel_to_fam)

svd.1000G <- bed_autoSVD(bed.1000G,
                         ind.col = info_snp$`_NUM_ID_.ss`, k = 30,
                         ncores = nb_cores())

library(ggplot2)
plot(svd.1000G) + scale_y_log10()
plot(svd.1000G, type = "scores")
plot(svd.1000G, type = "scores", scores = 1:30, cols = 5, coeff = 0.7)
plot(svd.1000G, type = "scores", scores = 17:18) + aes(color = is_rel)
plot(svd.1000G, type = "scores", scores = 21:22) + aes(color = is_rel)
plot(svd.1000G, type = "scores", scores = 23:24) + aes(color = is_rel)
plot(svd.1000G, type = "scores", scores = 25:26) + aes(color = is_rel)
plot(svd.1000G, type = "loadings", loadings = 1:10, coeff = 0.4)

dist <- robust::covRob(svd.1000G$u[, 15:30], estim = "pairwiseGK")$dist
plot(svd.1000G, type = "scores", scores = 17:18) +
  aes(color = log(dist)) +
  scale_colour_viridis_c()
plot(svd.1000G, type = "scores", scores = 21:22) +
  aes(color = log(dist)) +
  scale_colour_viridis_c()
plot(svd.1000G, type = "scores", scores = 23:24) +
  aes(color = log(dist)) +
  scale_colour_viridis_c()
plot(svd.1000G, type = "scores", scores = 25:26) +
  aes(color = log(dist)) +
  scale_colour_viridis_c()

which(log(dist) > 6)
hist(log(dist))

svd.1000G <- bed_autoSVD(bed.1000G, ind.row = which(log(dist) < 6),
                         ind.col = info_snp$`_NUM_ID_.ss`, k = 30,
                         ncores = nb_cores())

library(ggplot2)
plot(svd.1000G) + scale_y_log10()
plot(svd.1000G, type = "scores")
plot(svd.1000G, type = "scores", scores = 1:30, cols = 5, coeff = 0.7)


dist <- robust::covRob(svd.1000G$u[, 19:30], estim = "pairwiseGK")$dist
which(log(dist) > 6)
hist(log(dist))

svd.1000G <- bed_autoSVD(bed.1000G, ind.row = which(log(dist) < 6),
                         ind.col = info_snp$`_NUM_ID_.ss`, k = 30,
                         ncores = nb_cores())

library(ggplot2)
plot(svd.1000G) + scale_y_log10()
plot(svd.1000G, type = "scores")
plot(svd.1000G, type = "scores", scores = 1:30, cols = 5, coeff = 0.7)
plot(svd.1000G, type = "scores", scores = 17:18) + aes(color = is_rel)
plot(svd.1000G, type = "scores", scores = 21:22) + aes(color = is_rel)
plot(svd.1000G, type = "scores", scores = 23:24) + aes(color = is_rel)
plot(svd.1000G, type = "scores", scores = 25:26) + aes(color = is_rel)
plot(svd.1000G, type = "loadings", loadings = 1:10, coeff = 0.4)

dist <- robust::covRob(svd.1000G$u[, 15:30], estim = "pairwiseGK")$dist
plot(svd.1000G, type = "scores", scores = 17:18) +
  aes(color = log(dist)) +
  scale_colour_viridis_c()
plot(svd.1000G, type = "scores", scores = 21:22) +
  aes(color = log(dist)) +
  scale_colour_viridis_c()
plot(svd.1000G, type = "scores", scores = 23:24) +
  aes(color = log(dist)) +
  scale_colour_viridis_c()
plot(svd.1000G, type = "scores", scores = 25:26) +
  aes(color = log(dist)) +
  scale_colour_viridis_c()
