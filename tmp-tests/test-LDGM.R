files <- gtools::mixedsort(
  list.files("../datasets/EUR", full.names = TRUE))
length(files)  # 1361 -> number of blocks

snp_list <- gtools::mixedsort(
  list.files("../datasets/snplist", full.names = TRUE))


prec_info <- bigreadr::fread2(files[1])
snp_info <- bigreadr::fread2(snp_list[1])
prec <- Matrix::sparseMatrix(i = prec_info$V1, j = prec_info$V2, x = prec_info$V3,
                             dims = rep(nrow(snp_info), 2), index1 = FALSE,
                             symmetric = TRUE)
prec[1:5, 1:5]
hist(snp_info$EUR[Matrix::diag(prec) == 0], "FD")
hist(snp_info$EUR[Matrix::diag(prec) != 0], "FD")
