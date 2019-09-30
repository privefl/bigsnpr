################################################################################

#' Download 1000G
#'
#' Download 1000 genomes project (phase 3) data in PLINK bed/bim/fam format,
#' including 2493 (mostly unrelated) individuals
#' and ~1.4M SNPs in common with HapMap3.
#'
#' @param dir The directory where to put the downloaded files.
#' @param overwrite Whether to overwrite files when downloading and unzipping?
#'   Default is `FALSE`.
#' @param delete_zip Whether to delete zip after decompressing the file in it?
#'   Default is `TRUE`.
#'
#' @return The path of the downloaded bed file.
#'
#' @importFrom stats setNames
#'
#' @export
#'
download_1000G <- function(dir, overwrite = FALSE, delete_zip = TRUE) {

  files_unzipped <- file.path(dir, paste0("1000G_phase3_common_norel",
                                          c(".bed", ".bim", ".fam", ".fam2")))

  if (overwrite || !all(file.exists(files_unzipped))) {

    zip <- file.path(dir, "1000G_phase3_common_norel.zip")
    if (overwrite || !file.exists(zip)) {
      utils::download.file("https://ndownloader.figshare.com/files/17838962",
                           destfile = zip)
      if (delete_zip) on.exit(unlink(zip), add = TRUE)
    }
    utils::unzip(zip, exdir = dir)

  }

  files_unzipped[1]
}

################################################################################

#' Projecting PCA
#'
#' Projecting PCA of reference dataset to a target dataset.
#'
#' @param bed.new Object of type bed, which is the mapping of the bed file of
#'   the target data. Use `obj.bed <- bed(bedfile)` to get this object.
#' @param bed.ref Object of type bed, which is the mapping of the bed file of
#'   the reference data. Use `obj.bed <- bed(bedfile)` to get this object.
#' @param k Number of principal components to compute and project.
#' @param ind.row Rows to be used in the reference data. Default uses them all.
#' @inheritParams snp_match
#' @param match.min.prop Minimum proportion of variants in target data to be
#'   matched, otherwise stops with an error. Default is `50%`.
#' @param build.new Genome build of the target data. Default is `hg19`.
#' @param build.ref Genome build of the reference data. Default is `hg19`.
#' @inheritParams snp_modifyBuild
#' @param adjust String specifying the estimation method.
#'   Possible values are `"d.gsp"` (default), `"l.gsp"`, `"osp"`,
#'   or `"none"` (shrinkage-adjustement is skipped).
#' @inheritParams bed_autoSVD2
#' @inheritDotParams bed_autoSVD2 -obj.bed -ind.row -ind.col -k -verbose -ncores
#'
#' @return The projected `k` principal components.
#'   See `attr(<object>, "shrinkage")` to get the shrinkage values and
#'   `attr(<object>, "obj.svd.ref")` to get the partial decomposition of the
#'    reference data.
#' @export
#'
#' @seealso [hdpca::pc_adjust()]
#'
bed_projectPCA <- function(bed.new, bed.ref, k = 10,
                           ind.row = rows_along(bed.ref),
                           strand_flip = TRUE,
                           join_by_pos = TRUE,
                           match.min.prop = 0.5,
                           build.new = "hg19",
                           build.ref = "hg19",
                           liftOver = NULL,
                           ...,
                           adjust = c("d.gsp", "l.gsp", "osp", "none"),
                           verbose = TRUE,
                           ncores = 1) {

  adjust <- match.arg(adjust)

  # Verbose?
  printf2 <- function(...) if (verbose) printf(...)

  printf2("\n[Step 1/4] Matching variants of reference with target data..\n")

  map.ref <- setNames(bed.ref$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
  map.new <- setNames(bed.new$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
  if (join_by_pos && (build.new != build.ref)) {
    if (is.null(liftOver)) {
      stop2("Please provide path to liftOver executable.")
    } else if (substr(liftOver, 1, 1) != ".") {
      liftOver <- file.path(".", liftOver)
    }
    map.new <- snp_modifyBuild(map.new, liftOver = liftOver,
                               from = build.new, to = build.ref)
  }
  info_snp <- snp_match(cbind(map.ref, beta = 1), map.new,
                        strand_flip = strand_flip, join_by_pos = join_by_pos)
  if (nrow(info_snp) < (ncol(bed.new) * match.min.prop))
    stop2("Not enough variants have been matched.")

  printf2("\n[Step 2/4] Computing (auto) SVD of reference..\n")

  obj.svd <- bed_autoSVD2(bed.ref,
                          ind.row = ind.row,
                          ind.col = info_snp$`_NUM_ID_.ss`,
                          k = k,
                          ncores = ncores,
                          verbose = verbose,
                          ...)

  printf2("\n[Step 3/4] Projecting PC scores on new data..\n")

  keep <- match(attr(obj.svd, "subset.col"), info_snp$`_NUM_ID_.ss`)

  U2 <- big_apply(bed.new, function(X, ind, ind.col, center, scale, V) {
    ind.row <- rows_along(X)
    nb_nona <- bed_stats(X, ind.row, ind.col[ind])$nb_nona_row
    U <- apply(V[ind, , drop = FALSE], 2, function(v) {
      pMatVec4(X, ind.row, ind.col[ind], center[ind], scale[ind], v)
    })
    list(U, nb_nona)
  }, ind = seq_along(keep), ncores = ncores,
  ind.col = info_snp$`_NUM_ID_`[keep],
  center = (obj.svd$center - 1) * info_snp$beta[keep] + 1,
  scale = obj.svd$scale * info_snp$beta[keep],
  V = obj.svd$v)

  U3 <- U2[[1]][[1]]
  nb_nona <- U2[[1]][[2]]
  for (u2 in U2[-1]) {
    U3 <- U3 + u2[[1]]
    nb_nona <- nb_nona + u2[[2]]
  }
  U3 <- sweep(U3, 1, length(keep) / nb_nona, '*')

  printf2("\n[Step 4/4] Adjusting projected PC scores for shrinkage-bias..\n")

  if (adjust == "none") {
    printf2("Skipping adjustment.\n")
  } else {
    printf2("Computing all eigenvalues of reference..\n")
    K <- bed_tcrossprodSelf(bed.ref,
                            ind.row = attr(obj.svd, "subset.row"),
                            ind.col = attr(obj.svd, "subset.col"))
    eig.val <- eigen(K[], symmetric = TRUE, only.values = TRUE)$values
    nspike <- bigutilsr::pca_nspike(eig.val)
    printf2("Detected approximately %d 'distant spikes'.\n", nspike)
    U3 <- bigutilsr::pca_adjust(U3, eig.val, n.spikes = nspike,
                                p = length(attr(obj.svd, "subset.col")))
  }

  structure(U3, obj.svd.ref = obj.svd)
}

################################################################################
