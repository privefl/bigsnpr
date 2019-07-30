library(bigsnpr)

obj.bed <- bed.1000G <- bed("tmp-data/1000G_phase3_common_hapmap_norel.bed")
map.1000G <- bigreadr::fread2(
  bed.1000G$bimfile,
  col.names = c("chr", "rsid", "osef", "pos", "a1", "a0")
)

# bed_celiac <- "../paper2-PRS/backingfiles/celiacQC.bed"
bed_celiac <- "../Dubois2010_data/FinnuncorrNLITUK3hap550.bed"
info_snp <- snp_match(
  cbind(map.1000G[-3], beta = 1),
  setNames(bed(bed_celiac)$map[-3], c("chr", "rsid", "pos", "a1", "a0")),
  join_by_pos = FALSE
)
info_snp


{
  fam <- bed.1000G$fam
  ped <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
  fam2 <- dplyr::left_join(fam[c(2, 5)], ped[c(1:5, 7)],
                           by = c("sample.ID" = "Individual ID", "sex" = "Gender"))
  pop <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")
  fam2 <- dplyr::left_join(fam2, pop[1:3], by = c("Population" = "Population Code"))
  str(fam2)
}

ind <- which(fam2$`Super Population` == "EUR")
ind <- rows_along(fam2)


# AutoSVD algo
NCORES <- nb_cores()
nPC <- 20
ind.col <- info_snp$`_NUM_ID_.ss`
CHR <- infos.chr <- obj.bed$map$chromosome
POS <- infos.pos <- obj.bed$map$physical.pos

# Auto SVD
svd.1000G <- obj.svd <- bed_autoSVD2(obj.bed, ind.row = ind, ind.col = ind.col,
                                     k = nPC, ncores = NCORES)


library(ggplot2)
plot(svd.1000G)
plot(svd.1000G, type = "scores", scores = 1:20)
plot(svd.1000G, type = "scores") +
  aes(color = fam2$Population[attr(svd.1000G, "subset.row")]) +
  labs(color = "POP")
plot(svd.1000G, type = "scores", scores = 3:4) +
  aes(color = fam2$Population[attr(svd.1000G, "subset.row")]) +
  labs(color = "POP")

# snp_readBed(bed_celiac, "tmp-data/celiacQC2")
celiac <- snp_attach("tmp-data/celiacQC2.rds")
G <- celiac$genotypes

keep <- match(attr(svd.1000G, "subset.col"), info_snp$`_NUM_ID_.ss`)
options(bigstatsr.check.args = FALSE)
# why need to inverse beta?
U2 <- big_prodMat(G, svd.1000G$v,
                  ind.col = info_snp$`_NUM_ID_`[keep],
                  center = (svd.1000G$center - 1) * -info_snp$beta[keep] + 1,
                  scale = svd.1000G$scale * -info_snp$beta[keep])
dim(U2)  # 11950 x 20

obj.bed.celiac <- bed(bed_celiac)
U2 <- big_apply(obj.bed.celiac, function(X, ind, ind.row, ind.col, center, scale, V) {
  print(min(ind))
  nb_nona <- bigsnpr:::bed_stats(X, ind.row, ind.col[ind])$nb_nona_row
  U <- apply(V[ind, , drop = FALSE], 2, function(v) {
    bigsnpr:::pMatVec4(X, ind.row, ind.col[ind], center[ind], scale[ind], v)
  })
  list(U, nb_nona)
}, ind = seq_along(keep), ncores = nb_cores(),
ind.row = rows_along(obj.bed.celiac), ind.col = info_snp$`_NUM_ID_`[keep],
center = (svd.1000G$center - 1) * info_snp$beta[keep] + 1,
scale = svd.1000G$scale * info_snp$beta[keep],
V = svd.1000G$v)

U3 <- U2[[1]][[1]]
nona <- U2[[1]][[2]]
for (u2 in U2[-1]) {
  U3 <- U3 + u2[[1]]
  nona <- nona + u2[[2]]
}
U2 <- sweep(U3, 1, length(keep) / nona, '*')
dim(U2)

plot_grid(
  ggplot() +
    geom_point(aes(V1, V2), data = as.data.frame(predict(svd.1000G)), color = "red") +
    geom_point(aes(V1, V2), data = as.data.frame(U2), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V3, V4), data = as.data.frame(predict(svd.1000G)), color = "red") +
    geom_point(aes(V3, V4), data = as.data.frame(U2), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V5, V6), data = as.data.frame(predict(svd.1000G)), color = "red") +
    geom_point(aes(V5, V6), data = as.data.frame(U2), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V7, V8), data = as.data.frame(predict(svd.1000G)), color = "red") +
    geom_point(aes(V7, V8), data = as.data.frame(U2), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  scale = 0.9
)

# snp_readBed("tmp-data/1000G_phase3_common_hapmap_norel.bed")
snp_1000G <- snp_attach("tmp-data/1000G_phase3_common_hapmap_norel.rds")
G2 <- snp_1000G$genotypes
K <- big_tcrossprodSelf(G2, snp_scaleBinom(),
                        ind.row = attr(obj.svd, "subset.row"),
                        ind.col = attr(obj.svd, "subset.col"))
dim(K)
eigs <- eigen(K[], symmetric = TRUE, only.values = TRUE)
rbind(sqrt(head(eigs$values, length(obj.svd$d))), obj.svd$d)
# all.equal(attr(K, "center"), obj.svd$center)
# all.equal(attr(K, "scale"), obj.svd$scale)

train.eval <- eigs$values  # no sqrt()
p <- length(attr(obj.svd, "subset.col"))
n <- nrow(K)
testscore <- U2

score.adj.d2 <- hdpca::pc_adjust(train.eval, p, n, testscore,
                                 method = "d.gsp", n.spikes.max = nPC)
# score.adj.o2 <- hdpca::pc_adjust(train.eval,p,n,testscore,method="osp",n.spikes.max=10)
# score.adj.l2 <- hdpca::pc_adjust(train.eval,p,n,testscore,method="l.gsp",n.spikes.max=10)

pop <- fam2$`Super Population`[attr(svd.1000G, "subset.row")]
ggplot() +
  geom_point(aes(V1, V2, color = pop), data = as.data.frame(predict(svd.1000G))) +
  labs(color = "Pop 1000G") +
  geom_point(aes(V1, V2), data = as.data.frame(score.adj.d2), alpha = 0.1) +
  theme_bigstatsr()

ggplot() +
  geom_point(aes(V3, V4, color = pop), data = as.data.frame(predict(svd.1000G))) +
  labs(color = "Pop 1000G") +
  geom_point(aes(V3, V4), data = as.data.frame(score.adj.d2), alpha = 0.1) +
  theme_bigstatsr()

plot_grid(
  ggplot() +
    geom_point(aes(V1, V2), data = as.data.frame(predict(svd.1000G)), color = "red") +
    geom_point(aes(V1, V2), data = as.data.frame(score.adj.d2),
               color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V3, V4), data = as.data.frame(predict(svd.1000G)), color = "red") +
    geom_point(aes(V3, V4), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V5, V6), data = as.data.frame(predict(svd.1000G)), color = "red") +
    geom_point(aes(V5, V6), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V7, V8), data = as.data.frame(predict(svd.1000G)), color = "red") +
    geom_point(aes(V7, V8), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  scale = 0.9
)

plot_grid(
  ggplot() +
    geom_point(aes(V1, V2), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V1, V2), data = as.data.frame(predict(svd.1000G)), color = "red") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V3, V4), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V3, V4), data = as.data.frame(predict(svd.1000G)), color = "red") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V5, V6), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V5, V6), data = as.data.frame(predict(svd.1000G)), color = "red") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V7, V8), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V7, V8), data = as.data.frame(predict(svd.1000G)), color = "red") +
    theme_bigstatsr(),

  scale = 0.9
)

plot_grid(
  ggplot() +
    geom_point(aes(V9, V10), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V9, V10), data = as.data.frame(predict(svd.1000G)), color = "red") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V11, V12), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V11, V12), data = as.data.frame(predict(svd.1000G)), color = "red") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V13, V14), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V13, V14), data = as.data.frame(predict(svd.1000G)), color = "red") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V15, V16), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V15, V16), data = as.data.frame(predict(svd.1000G)), color = "red") +
    theme_bigstatsr(),

  scale = 0.9
)

# obj.dist <- robust::covRob(predict(svd.1000G)[pop == "EUR", c(1:4)],
#                            estim = "pairwiseGK")
# dist <- stats::mahalanobis(score.adj.d2[, c(1:4)], obj.dist$center, obj.dist$cov)
# dist <- dbscan::lof(score.adj.d2[, 1:4], k = 4)
#
# # svd.celiac <- snp_autoSVD(G, celiac$map$chromosome, celiac$map$physical.pos,
# #                           thr.r2 = 0.1, ncores = NCORES)
# # saveRDS(svd.celiac, "tmp-data/svd_celiac.rds")
# svd.celiac <- readRDS("tmp-data/svd_celiac.rds")
# plot(svd.celiac)
#
# plot(svd.celiac, type = "scores", scores = 5:6) +
#   aes(color = log(dist)) + scale_colour_viridis_c()
