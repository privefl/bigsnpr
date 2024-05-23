################################################################################

context("CORR_EXT_THR")

################################################################################

################################################################################

test_that("snp_cor_extendedThr (no subset) works", {

  bigsnp <- snp_attachExtdata()
  G <- bigsnp$genotypes
  POS <- bigsnp$map$physical.pos
  corr0 <- big_cor(G)

  replicate(5, {

    THR <- runif(1, 0.01, 0.5)
    SIZE <- round(runif(1, 100, 2000))
    test <- snp_cor_extendedThr(G, thr_r2 = THR, infos.pos = POS, size = SIZE, ncores = NCORES)

    corr_T <- as(Matrix::drop0(test, tol = sqrt(THR)), "generalMatrix")
    list_keep <- bigsnpr:::find_indirect_corr(corr_T@p, corr_T@i, ncores = NCORES)
    overlap <- Matrix::which((corr_T %*% corr_T) != 0, arr.ind = TRUE)
    check_list <- split(overlap[, "row"], factor(cols_along(G))[overlap[, "col"]])

    expect_equal(list_keep, unname(check_list))
    expect_equal(test[overlap], corr0[overlap])

    NULL
  })
})

################################################################################

test_that("snp_cor_extendedThr (no subset) works", {

  bigsnp <- snp_attachExtdata()
  G <- bigsnp$genotypes
  POS <- bigsnp$map$physical.pos

  replicate(5, {

    THR <- runif(1, 0.01, 0.5)
    SIZE <- round(runif(1, 100, 2000))
    ind.row <- sample(nrow(G), 500)
    ind.col <- sort(sample(ncol(G), 1000))

    G2 <- big_copy(G, ind.row = ind.row, ind.col = ind.col)
    G2[1, ] <- 3  # setting as missing, will not be used
    test0 <- snp_cor_extendedThr(G2, thr_r2 = THR, infos.pos = POS[ind.col],
                                 size = SIZE, ncores = NCORES)

    test <- snp_cor_extendedThr(G, thr_r2 = THR, infos.pos = POS[ind.col],
                                size = SIZE, ncores = NCORES,
                                ind.row = ind.row[-1], ind.col = ind.col)

    expect_equal(test, test0)

    NULL
  })
})

################################################################################
