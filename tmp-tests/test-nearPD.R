corr <- runonce::save_run(
  snp_cor(chr22$genotypes, infos.pos = POS2, size = 3 / 1000, ncores = 6),
  file = "tmp-data/corr_chr22.rds"
)

corr_sub <- corr[1:10, 1:10]
round(100 * corr_sub, 2)

options(max.print = 330)
round(100 * Matrix::nearPD(corr_sub, base.matrix = TRUE)$mat, 2)

corr_sub_thr <- Matrix::drop0(corr_sub, 0.1)  # still PD
k <- 6
(corr_sub_thr <- Matrix::band(corr_sub, -k, k))
eigen(corr_sub_thr, symmetric = TRUE, only.values = TRUE)$val
round(100 * corr_sub_thr, 2)
round(100 * Matrix::nearPD(corr_sub_thr, base.matrix = TRUE,
                           corr = TRUE, doSym = TRUE)$mat, 2)
# can get larger values..
round(100 * corr_sub, 2)
