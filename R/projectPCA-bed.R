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

part_prod <- function(X, ind, ind.row, ind.col, center, scale, V, XV, X_norm) {

  res <- prod_and_rowSumsSq(
    obj_bed = X,
    ind_row = ind.row,
    ind_col = ind.col[ind],
    center  = center[ind],
    scale   = scale[ind],
    V       = V[ind, , drop = FALSE]
  )

  big_increment(XV,     res[[1]], use_lock = TRUE)
  big_increment(X_norm, res[[2]], use_lock = TRUE)
}

################################################################################

#' Projecting PCA
#'
#' Computing and projecting PCA of reference dataset to a target dataset.
#'
#' @param obj.bed.new Object of type bed, which is the mapping of the bed file of
#'   the target data. Use `obj.bed <- bed(bedfile)` to get this object.
#' @param obj.bed.ref Object of type bed, which is the mapping of the bed file of
#'   the reference data. Use `obj.bed <- bed(bedfile)` to get this object.
#' @param k Number of principal components to compute and project.
#' @param ind.row.new Rows to be used in the target data. Default uses them all.
#' @param ind.row.ref Rows to be used in the reference data.
#'   Default uses them all.
#' @param ind.col.ref Columns to be potentially used in the reference data.
#'   Default uses all the ones in common with target data.
#' @inheritParams snp_match
#' @param match.min.prop Minimum proportion of variants in target data to be
#'   matched, otherwise stops with an error. Default is `50%`.
#' @param build.new Genome build of the target data. Default is `hg19`.
#' @param build.ref Genome build of the reference data. Default is `hg19`.
#' @inheritParams snp_modifyBuild
#' @inheritParams bed_autoSVD2
#' @inheritDotParams bed_autoSVD2 -obj.bed -ind.row -ind.col -k -verbose -ncores
#'
#' @return A list of 3 elements:
#'   - `$obj.svd.ref`: big_SVD object computed from reference data.
#'   - `$simple_proj`: simple projection of new data into space of reference PCA.
#'   - `$OADP_proj`: Online Augmentation, Decomposition, and Procrustes (OADP)
#'     projection of new data into space of reference PCA.
#'
#' @export
#'
bed_projectPCA <- function(obj.bed.ref, obj.bed.new, k = 10,
                           ind.row.new = rows_along(obj.bed.new),
                           ind.row.ref = rows_along(obj.bed.ref),
                           ind.col.ref = cols_along(obj.bed.ref),
                           strand_flip = TRUE,
                           join_by_pos = TRUE,
                           match.min.prop = 0.5,
                           build.new = "hg19",
                           build.ref = "hg19",
                           liftOver = NULL,
                           ...,
                           verbose = TRUE,
                           ncores = 1) {

  # Verbose?
  printf2 <- function(...) if (verbose) printf(...)

  printf2("\n[Step 1/3] Matching variants of reference with target data..\n")

  map.ref <- setNames(obj.bed.ref$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
  map.new <- setNames(obj.bed.new$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
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
  if (nrow(info_snp) < (ncol(obj.bed.new) * match.min.prop))
    stop2("Not enough variants have been matched.")

  printf2("\n[Step 2/3] Computing (auto) SVD of reference..\n")

  obj.svd <- bed_autoSVD2(
    obj.bed.ref,
    ind.row = ind.row.ref,
    ind.col = intersect(ind.col.ref, info_snp$`_NUM_ID_.ss`),
    k       = k,
    ncores  = ncores,
    verbose = verbose,
    ...
  )

  printf2("\n[Step 3/3] Projecting PC scores on new data..\n")

  keep   <- match(attr(obj.svd, "subset"), info_snp$`_NUM_ID_.ss`)
  X_norm <- FBM(length(ind.row.new), 1, init = 0)
  XV     <- FBM(length(ind.row.new), k, init = 0)

  big_parallelize(
    obj.bed.new,
    p.FUN = part_prod,
    ind = seq_along(keep),
    ncores = ncores,
    ind.row = ind.row.new,
    ind.col = info_snp$`_NUM_ID_`[keep],
    center = (obj.svd$center - 1) * info_snp$beta[keep] + 1,
    scale = obj.svd$scale * info_snp$beta[keep],
    V = obj.svd$v,
    XV = XV,
    X_norm = X_norm
  )

  XV <- XV[]

  list(
    obj.svd.ref = obj.svd,
    simple_proj = XV,
    OADP_proj   = bigutilsr::pca_OADP_proj2(XV, X_norm[], obj.svd$d)
  )
}

################################################################################

#' Projecting PCA
#'
#' Projecting PCA using individuals from one dataset
#' to other individuals from the same dataset.
#'
#' @param obj.svd List with `v`, `d`, `center` and `scale`. Typically the an
#'   object of type "big_SVD".
#' @param obj.bed Object of type bed, which is the mapping of the bed file of
#'   the data containing both the individuals that were used to compute the PCA
#'   and the other individuals to be projected.
#' @param ind.row Rows (individuals) to be projected.
#' @param ind.col Columns that were used for computing PCA. If [bed_autoSVD2] was
#'   used, then `attr(obj.svd, "subset")` is automatically used by default.
#'   Otherwise (e.g. if [bed_randomSVD] was used), you have to pass `ind.col`.
#' @inheritParams bigstatsr::big_parallelize
#'
#' @inherit bed_projectPCA return
#'
#' @export
#'
bed_projectSelfPCA <- function(obj.svd, obj.bed, ind.row,
                               ind.col = attr(obj.svd, "subset"),
                               ncores = 1) {

  assert_lengths(rows_along(obj.svd$v), ind.col)

  X_norm <- FBM(length(ind.row), 1,               init = 0)
  XV     <- FBM(length(ind.row), ncol(obj.svd$v), init = 0)

  big_parallelize(
    obj.bed,
    p.FUN = part_prod,
    ind = seq_along(ind.col),
    ncores = ncores,
    ind.row = ind.row,
    ind.col = ind.col,
    center = obj.svd$center,
    scale = obj.svd$scale,
    V = obj.svd$v,
    XV = XV,
    X_norm = X_norm
  )

  XV <- XV[]

  list(
    obj.svd.ref = obj.svd,
    simple_proj = XV,
    OADP_proj   = bigutilsr::pca_OADP_proj2(XV, X_norm[], obj.svd$d)
  )
}

################################################################################
