################################################################################

context("PCADAPT") # Need to change tolerances

################################################################################

bigsnp <- snp_attachExtdata()
expect_s3_class(bigsnp, "bigSNP")
G <- bigsnp$genotypes
expect_s4_class(G, "FBM.code256")

################################################################################

test_that("Same as pcadapt", {

  # skip_on_cran()

  ################################################################################

  tmpfile <- tempfile(fileext = ".pcadapt")
  write.table(t(G[]), tmpfile, quote = FALSE, sep = " ",
              row.names = FALSE, col.names = FALSE)

  expect_output(bed <- pcadapt::read.pcadapt(tmpfile, type = "pcadapt"))
  obj.pcadapt <- pcadapt::pcadapt(bed, K = 10, min.maf = 0)

  ################################################################################

  obj.svd <- big_SVD(G, snp_scaleBinom())
  test <- bigsnpr:::multLinReg(G, rows_along(G), cols_along(G), obj.svd$u)

  expect_equal(obj.svd$center / 2, obj.pcadapt$af)
  expect_equal(obj.svd$d, obj.pcadapt$singular.values * sqrt((nrow(G) - 1) * ncol(G)))
  expect_equal(abs(obj.svd$u), abs(obj.pcadapt$scores),   tolerance = 1e-4)
  expect_equal(abs(obj.svd$v), abs(obj.pcadapt$loadings), tolerance = 1e-2)
  expect_equal(test,           obj.pcadapt$zscores,       tolerance = 0.2)
  # expect_equal(abs(cov(obj.svd$u, obj.pcadapt$scores)), diag(10) / (nrow(G) - 1))
  # expect_equal(abs(cor(obj.svd$v, obj.pcadapt$loadings)), diag(10),
  #              tolerance = 1e-2, check.attributes = FALSE)
  # expect_equal(abs(cor(test, obj.pcadapt$zscores)), diag(10),
  #              tolerance = 1e-2, check.attributes = FALSE)

  ################################################################################

  obj.gwas <- snp_pcadapt(G, obj.svd$u, ncores = sample(1:2, 1))
  expect_equal(bigsnpr:::getLambdaGC(obj.gwas), 1)  ## always GC
  expect_equal(get("lamGC", envir = environment(attr(obj.gwas, "transfo"))),
               obj.pcadapt$gif, tolerance = 1e-2)
  expect_equal(obj.gwas$score, as.numeric(obj.pcadapt$stat), tolerance = 1e-1)
  expect_equal(predict(obj.gwas, log10 = FALSE), obj.pcadapt$pvalues,
               tolerance = 1e-1)

  ################################################################################

  expect_s3_class(snp_qq(obj.gwas, lambdaGC = FALSE), "ggplot")
  p <- snp_qq(snp_gc(obj.gwas))
  expect_equal(as.character(p$labels$subtitle), "lambda[GC] == 1")
  expect_s3_class(p, "ggplot")
  p2 <- snp_manhattan(obj.gwas, infos.chr = bigsnp$map$chromosome,
                      infos.pos = bigsnp$map$physical.pos,
                      npoints = 2000)
  expect_s3_class(p2, "ggplot")

  ################################################################################

  # K = 1
  obj.pcadapt <- pcadapt::pcadapt(bed, K = 1, min.maf = 0)
  obj.gwas    <- snp_pcadapt(G, obj.svd$u[, 1], ncores = sample(1:2, 1))

  expect_s3_class(snp_qq(obj.gwas), "ggplot")

  obj.bed <- bed(snp_writeBed(bigsnp, bedfile = tempfile(fileext = ".bed")))
  obj.gwas2 <- bed_pcadapt(obj.bed, obj.svd$u[, 1], ncores = sample(1:2, 1))
  expect_equal(obj.gwas2, obj.gwas)

  expect_equal(bigsnpr:::getLambdaGC(obj.gwas), 1)  ## always GC
  expect_equal(get("lamGC", envir = environment(attr(obj.gwas, "transfo"))),
               obj.pcadapt$gif, tolerance = 1e-2)
  expect_equal(obj.gwas$score, as.numeric(obj.pcadapt$stat), tolerance = 1e-2)
  expect_equal(predict(obj.gwas, log10 = FALSE), obj.pcadapt$pvalues,
               tolerance = 1e-2)
  expect_equal(obj.pcadapt$scores[, 1],   obj.svd$u[, 1], tolerance = 1e-8)
  expect_equal(obj.pcadapt$loadings[, 1], obj.svd$v[, 1], tolerance = 1e-8)

  ################################################################################

})
