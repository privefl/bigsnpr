library(bigsnpr)

bedfile.new <- "../Dubois2010_data/FinnuncorrNLITUK3hap550.bed"
bed.new <- bed(bedfile.new)
bedfile.ref <- "tmp-data/1000G_phase3_common_hapmap_norel.bed"
bed.ref <- bed(bedfile.ref)
ind.row <- rows_along(bed.ref)
strand_flip = TRUE
join_by_pos = FALSE  # TRUE by default
ncores <- nb_cores()
nPC <- 20

adjust <- c("d.gsp", "l.gsp", "osp", "none")
adjust <- adjust[1] # match.arg(adjust)

verbose <- TRUE
# Verbose?
printf2 <- function(...) if (verbose) bigsnpr:::printf(...)

printf2("[Step 1/5] Matching variants with reference..\n")

ref.map <- setNames(bed.ref$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
map     <- setNames(bed.new$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
info_snp <- snp_match(cbind(ref.map, beta = 1), map,
                      strand_flip = strand_flip,
                      join_by_pos = join_by_pos)


# Auto SVD -> TODO: add '...'
printf2("[Step 2/5] Computing (auto) SVD of reference..\n")

obj.svd <- bed_autoSVD2(bed.ref,
                        ind.row = ind.row,
                        ind.col = info_snp$`_NUM_ID_.ss`,
                        k = nPC,
                        ncores = ncores)


printf2("[Step 3/5] Projecting PC scores on new data..\n")

keep <- match(attr(obj.svd, "subset.col"), info_snp$`_NUM_ID_.ss`)

U2 <- big_apply(bed.new, function(X, ind, ind.col, center, scale, V) {
  ind.row <- rows_along(X)
  nb_nona <- bigsnpr:::bed_stats(X, ind.row, ind.col[ind])$nb_nona_row
  U <- apply(V[ind, , drop = FALSE], 2, function(v) {
    bigsnpr:::pMatVec4(X, ind.row, ind.col[ind], center[ind], scale[ind], v)
  })
  list(U, nb_nona)
}, ind = seq_along(keep), ncores = ncores,
ind.col = info_snp$`_NUM_ID_`[keep],
center = (obj.svd$center - 1) * info_snp$beta[keep] + 1,
scale = obj.svd$scale * info_snp$beta[keep],
V = obj.svd$v)

U3 <- U2[[1]][[1]]
nona <- U2[[1]][[2]]
for (u2 in U2[-1]) {
  U3 <- U3 + u2[[1]]
  nona <- nona + u2[[2]]
}
U3 <- sweep(U3, 1, length(keep) / nona, '*')
# dim(U3)

library(ggplot2)
plot_grid(
  ggplot() +
    geom_point(aes(V1, V2), data = as.data.frame(predict(obj.svd)), color = "red") +
    geom_point(aes(V1, V2), data = as.data.frame(U3), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V3, V4), data = as.data.frame(predict(obj.svd)), color = "red") +
    geom_point(aes(V3, V4), data = as.data.frame(U3), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V5, V6), data = as.data.frame(predict(obj.svd)), color = "red") +
    geom_point(aes(V5, V6), data = as.data.frame(U3), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V7, V8), data = as.data.frame(predict(obj.svd)), color = "red") +
    geom_point(aes(V7, V8), data = as.data.frame(U3), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  scale = 0.9
)


printf2("[Step 4/5] Adjusting projected PC scores..\n")

K <- bed_tcrossprodSelf(bed.ref,
                        ind.row = attr(obj.svd, "subset.row"),
                        ind.col = attr(obj.svd, "subset.col"))
dim(K)
eig.val <- eigen(K[], symmetric = TRUE, only.values = TRUE)$values
rbind(sqrt(head(eig.val, length(obj.svd$d))), obj.svd$d)

p <- length(attr(obj.svd, "subset.col"))
n <- nrow(K)
testscore <- U3


system.time(test <- hdpca::select.nspike(eig.val, p, n, n.spikes.max = 50))
test
S <- eig.val / p
hist(S2 <- log(head(eig.val, -1)), breaks = 50, freq = FALSE)
curve(dnorm(x, median(S2), mad(S2)), col = "red", add = TRUE)
print(q <- bigsnpr:::tukey_MC_up(S2, coef = NULL))
abline(v = q, col = "red")
sum(S2 > q, na.rm = TRUE)
hist(S2[S2 > 12])

hist(pval <- pnorm(S2, median(S2), mad(S2), lower.tail = FALSE))
plot(head(pval, 80))

hist(S3 <- diff(diff(S2)), breaks = 50)
which(S3 > bigsnpr:::tukey_MC_up(S3))

system.time(
  score.adj.d2 <- hdpca::pc_adjust(eig.val, p, n, testscore,
                                   method = "d.gsp", n.spikes.max = nPC)
)

system.time(
  score.adj.o2 <- hdpca::pc_adjust(eig.val,p,n,testscore,
                                   method="osp",n.spikes.max=nPC)
)
all.equal(score.adj.o2, score.adj.d2)

system.time(
  score.adj.l2 <- hdpca::pc_adjust(eig.val,p,n,testscore,
                                   method="l.gsp",n.spikes.max=nPC)
)
all.equal(score.adj.l2, score.adj.o2)



{
  fam <- bed.ref$fam
  ped <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
  fam2 <- dplyr::left_join(fam[c(2, 5)], ped[c(1:5, 7)],
                           by = c("sample.ID" = "Individual ID", "sex" = "Gender"))
  pop <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")
  fam2 <- dplyr::left_join(fam2, pop[1:3], by = c("Population" = "Population Code"))
  str(fam2)
}

pop <- fam2$`Super Population`[attr(obj.svd, "subset.row")]
ggplot() +
  geom_point(aes(V1, V2, color = pop), data = as.data.frame(predict(obj.svd))) +
  labs(color = "Pop 1000G") +
  geom_point(aes(V1, V2), data = as.data.frame(score.adj.d2), alpha = 0.1) +
  theme_bigstatsr()

ggplot() +
  geom_point(aes(V3, V4, color = pop), data = as.data.frame(predict(obj.svd))) +
  labs(color = "Pop 1000G") +
  geom_point(aes(V3, V4), data = as.data.frame(score.adj.d2), alpha = 0.1) +
  theme_bigstatsr()

plot_grid(
  ggplot() +
    geom_point(aes(V1, V2), data = as.data.frame(predict(obj.svd)), color = "red") +
    geom_point(aes(V1, V2), data = as.data.frame(score.adj.d2),
               color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V3, V4), data = as.data.frame(predict(obj.svd)), color = "red") +
    geom_point(aes(V3, V4), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V5, V6), data = as.data.frame(predict(obj.svd)), color = "red") +
    geom_point(aes(V5, V6), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V7, V8), data = as.data.frame(predict(obj.svd)), color = "red") +
    geom_point(aes(V7, V8), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    theme_bigstatsr(),

  scale = 0.9
)

plot_grid(
  ggplot() +
    geom_point(aes(V1, V2), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V1, V2), data = as.data.frame(predict(obj.svd)), color = "red") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V3, V4), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V3, V4), data = as.data.frame(predict(obj.svd)), color = "red") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V5, V6), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V5, V6), data = as.data.frame(predict(obj.svd)), color = "red") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V7, V8), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V7, V8), data = as.data.frame(predict(obj.svd)), color = "red") +
    theme_bigstatsr(),

  scale = 0.9
)

plot_grid(
  ggplot() +
    geom_point(aes(V9, V10), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V9, V10), data = as.data.frame(predict(obj.svd)), color = "red") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V11, V12), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V11, V12), data = as.data.frame(predict(obj.svd)), color = "red") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V13, V14), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V13, V14), data = as.data.frame(predict(obj.svd)), color = "red") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V15, V16), data = as.data.frame(score.adj.d2), color = "blue", alpha = 0.3) +
    geom_point(aes(V15, V16), data = as.data.frame(predict(obj.svd)), color = "red") +
    theme_bigstatsr(),

  scale = 0.9
)



