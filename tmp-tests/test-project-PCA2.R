library(bigsnpr)

bedfile.new <- "../Dubois2010_data/FinnuncorrNLITUK3hap550.bed"; build <- "hg19"
bedfile.new <- "../POPRES_data/POPRES_allchr.bed"; build <- "hg18"
bed.new <- bed(bedfile.new)
bedfile.ref <- download_1000G("tmp-data")
bed.ref <- bed(bedfile.ref)
ind.row <- rows_along(bed.ref)
strand_flip = TRUE
join_by_pos = TRUE
ncores <- nb_cores()
nPC <- 20

adjust <- c("d.gsp", "l.gsp", "osp", "none")
adjust <- adjust[1] # match.arg(adjust)

verbose <- TRUE
# Verbose?
printf2 <- function(...) if (verbose) bigsnpr:::printf(...)

printf2("[Step 1/5] Matching variants with reference..\n")

ref.map <- setNames(bed.ref$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
if (build != "hg19") ref.map <- snp_modifyBuild(ref.map, "./tmp-data/liftOver",
                                                from = "hg19", to = build)
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
  ind.row <- bigstatsr::rows_along(X)
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
# 587 sec
test
test$n.spikes # 21
system.time(test2 <- bigutilsr::pca_nspike(eig.val))
# 2 sec
test2 # 21

system.time(
  score.adj.d2 <- bigutilsr::pca_adjust(testscore, eig.val, p, n.spikes = 8)
)
attr(score.adj.d2, "shrinkage")
# 0.9960032 0.9897782 0.9628025 0.9519429 0.7709360 0.7649171 0.7344308 0.7008195
# 0.5068197 0.4748159 0.4146323 0.3976913 0.3599203 0.3287099 0.3199303
# 0.9959816 0.9897231 0.9626026 0.9516850 0.7696716 0.7636153 0.7329297 0.6990734

fam2 <- bigreadr::fread2(sub_bed(bedfile.ref, ".fam2"))[attr(obj.svd, "subset.row"), ]

pop <- fam2$`Super Population`
ggplot() +
  geom_point(aes(V1, V2, color = pop), data = as.data.frame(predict(obj.svd))) +
  labs(color = "Pop 1000G") +
  geom_point(aes(V1, V2), data = as.data.frame(score.adj.d2), alpha = 0.3) +
  theme_bigstatsr()

ggplot() +
  geom_point(aes(V3, V4, color = pop), data = as.data.frame(predict(obj.svd))) +
  labs(color = "Pop 1000G") +
  geom_point(aes(V3, V4), data = as.data.frame(score.adj.d2), alpha = 0.3) +
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

plot_grid(

  ggplot() +
    geom_point(aes(V1, V2), data = as.data.frame(score.adj.d2)) +
    geom_point(aes(V1, V2, color = pop), alpha = 0.4,
               data = as.data.frame(predict(obj.svd))) +
    labs(color = "Pop 1000G") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V3, V4), data = as.data.frame(score.adj.d2)) +
    geom_point(aes(V3, V4, color = pop), alpha = 0.4,
               data = as.data.frame(predict(obj.svd))) +
    labs(color = "Pop 1000G") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V5, V6), data = as.data.frame(score.adj.d2)) +
    geom_point(aes(V5, V6, color = pop), alpha = 0.4,
               data = as.data.frame(predict(obj.svd))) +
    labs(color = "Pop 1000G") +
    theme_bigstatsr(),

  ggplot() +
    geom_point(aes(V7, V8), data = as.data.frame(score.adj.d2)) +
    geom_point(aes(V7, V8, color = pop), alpha = 0.4,
               data = as.data.frame(predict(obj.svd))) +
    labs(color = "Pop 1000G") +
    theme_bigstatsr()

)

hist(S <- bigutilsr::covRob(score.adj.d2, estim = "pairwiseGK")$dist, breaks = 50)
hist(S2 <- log(S), breaks = length(S) / 10)
abline(v = (q <- bigutilsr::tukey_mc_up(S2)), col = "red")
hist(S, breaks = length(S) / 10)
hist(log(S2), breaks = length(S) / 10)
bigutilsr::pca_nspike(log(S2))
hist(S[S < 60], breaks = length(S) / 10)
bigutilsr::pca_nspike(S)

