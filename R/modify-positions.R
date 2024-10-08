################################################################################

#' Modify genome build
#'
#' Modify the physical position information of a data frame
#' when converting genome build using executable *liftOver*.
#'
#' @param info_snp A data frame with columns "chr" and "pos".
#' @param liftOver Path to liftOver executable. Binaries can be downloaded at
#'   \url{https://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/liftOver} for Mac
#'   and at \url{https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver}
#'   for Linux.
#' @param from Genome build to convert from. Default is `hg18`.
#' @param to Genome build to convert to. Default is `hg19`.
#' @param check_reverse Whether to discard positions for which we cannot go back
#'   to initial values by doing 'from -> to -> from'. Default is `TRUE`.
#' @param local_chain Local chain file (e.g. `hg18ToHg19.over.chain.gz`) to use
#'   instead of downloading one from parameters `from` and `to` (the default).
#'   You can download one such file from e.g.
#'   \url{https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/}.
#'   Provide a vector of two when using `check_reverse`.
#' @param base_url From where to download the chain files. Default is
#'   `"https://hgdownload.soe.ucsc.edu/goldenPath/"`. You can also try
#'   replacing `https` by `http`, and/or `soe` by `cse`.
#'
#' @references
#' Hinrichs, Angela S., et al. "The UCSC genome browser database: update 2006."
#' Nucleic acids research 34.suppl_1 (2006): D590-D598.
#'
#' @return Input data frame `info_snp` with column "pos" in the new build.
#' @export
#'
snp_modifyBuild <- function(info_snp, liftOver,
                            from = "hg18", to = "hg19",
                            check_reverse = TRUE,
                            local_chain = NULL,
                            base_url = "https://hgdownload.soe.ucsc.edu/goldenPath/") {

  if (!all(c("chr", "pos") %in% names(info_snp)))
    stop2("Expecting variables 'chr' and 'pos' in input 'info_snp'.")

  # Make sure liftOver is executable
  liftOver <- make_executable(normalizePath(liftOver))

  # Need BED UCSC file for liftOver
  info_BED <- with(info_snp, data.frame(
    # sub("^0", "", c("01", 1, 22, "X")) -> "1"  "1"  "22" "X"
    chrom = paste0("chr", sub("^0", "", chr)),
    start = pos - 1L, end = pos,
    id = seq_along(pos)))

  BED <- tempfile(fileext = ".BED")
  bigreadr::fwrite2(stats::na.omit(info_BED),
                    BED, col.names = FALSE, sep = " ", scipen = 50)

  # Need chain file
  if (is.null(local_chain)) {
    url <- paste0(base_url, from, "/liftOver/",
                  from, "To", tools::toTitleCase(to), ".over.chain.gz")
    chain <- tempfile(fileext = ".over.chain.gz")
    utils::download.file(url, destfile = chain, quiet = TRUE)
  } else {
    chain <- local_chain[[1]]
    assert_exist(chain)
  }

  # Run liftOver (usage: liftOver oldFile map.chain newFile unMapped)
  lifted <- tempfile(fileext = ".BED")
  system2(liftOver, c(BED, chain, lifted, tempfile(fileext = ".txt")))

  # Read the ones lifter + some QC
  new_pos <- bigreadr::fread2(lifted, nThread = 1)
  is_bad <- vctrs::vec_duplicate_detect(new_pos$V4) |
    (new_pos$V1 != info_BED$chrom[new_pos$V4])
  new_pos <- new_pos[which(!is_bad), ]

  pos0 <- info_snp$pos
  info_snp$pos <- NA_integer_
  info_snp$pos[new_pos$V4] <- new_pos$V3

  if (check_reverse) {
    pos2 <- suppressMessages(
      Recall(info_snp, liftOver, from = to, to = from, check_reverse = FALSE,
             local_chain = local_chain[[2]], base_url = base_url)$pos
    )
    info_snp$pos[pos2 != pos0] <- NA_integer_
  }

  message2("%d variants have not been mapped.", sum(is.na(info_snp$pos)))

  info_snp
}

################################################################################

#' Interpolate to genetic positions
#'
#' Use genetic maps available at
#' \url{https://github.com/joepickrell/1000-genomes-genetic-maps/}
#' to interpolate physical positions (in bp) to genetic positions (in cM).
#'
#' @inheritParams bigsnpr-package
#' @param dir Directory where to download and decompress files.
#'   Default is `tempdir()`. Directly use *uncompressed* files there if already
#'   present. You can use [R.utils::gunzip()] to uncompress local files.
#' @param rsid If providing rsIDs, the matching is performed using those
#'   (instead of positions) and variants not matched are interpolated using
#'   spline interpolation of variants that have been matched.
#' @param type Whether to use the genetic maps interpolated from "OMNI"
#'   (the default), or from "hapmap".
#'
#' @return The new vector of genetic positions.
#' @export
#'
snp_asGeneticPos <- function(infos.chr, infos.pos, dir = tempdir(), ncores = 1,
                             rsid = NULL, type = c("OMNI", "hapmap")) {

  type <- match.arg(type)
  path <- c(OMNI   = "interpolated_OMNI",
            hapmap = "interpolated_from_hapmap")[type]

  assert_package("R.utils")
  assert_lengths(infos.chr, infos.pos)
  if (!is.null(rsid)) assert_lengths(rsid, infos.pos)

  snp_split(infos.chr, function(ind.chr, pos, dir, rsid) {

    chr <- attr(ind.chr, "chr")
    basename <- paste0("chr", chr, `if`(type == "OMNI", ".OMNI", ""),
                       ".interpolated_genetic_map")
    mapfile <- file.path(dir, basename)
    if (!file.exists(mapfile)) {
      url <- paste0("https://github.com/joepickrell/1000-genomes-genetic-maps/",
                    "raw/master/", path, "/", basename, ".gz")
      gzfile <- paste0(mapfile, ".gz")
      utils::download.file(url, destfile = gzfile, quiet = TRUE)
      R.utils::gunzip(gzfile)
    }
    map.chr <- bigreadr::fread2(mapfile, showProgress = FALSE, nThread = 1)

    if (is.null(rsid)) {
      ind <- bigutilsr::knn_parallel(as.matrix(map.chr$V2), as.matrix(pos[ind.chr]),
                                     k = 1, ncores = 1)$nn.idx
      new_pos <- map.chr$V3[ind]
    } else {
      ind <- match(rsid[ind.chr], map.chr$V1)
      new_pos <- map.chr$V3[ind]

      indNA <- which(is.na(ind))
      if (length(indNA) > 0) {
        pos.chr <- pos[ind.chr]
        new_pos[indNA] <- suppressWarnings(
          stats::spline(pos.chr, new_pos, xout = pos.chr[indNA], method = "hyman")$y)
      }
    }

    new_pos

  }, combine = "c", pos = infos.pos, dir = dir, rsid = rsid, ncores = ncores)
}

################################################################################

#' Download a genetic map
#'
#' @param type Which genetic map to download.
#' @param dir Directory where to download and decompress files.
#' @inheritParams bigsnpr-package
#'
#' @return A data frame with 3 columns: `chr`, `pos`, and `pos_cM`.
#' @export
#'
#' @details
#' The hg19 genetic maps are downloaded from
#' \url{https://github.com/joepickrell/1000-genomes-genetic-maps/}
#' while the hg38 one is downloaded from
#' `https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/`.
#'
#' @rdname snp_asGeneticPos2
#'
download_genetic_map <- function(type = c("hg19_OMNI", "hg19_hapmap", "hg38_price"),
                                 dir, ncores = 1) {

  assert_package("R.utils")

  type <- match.arg(type)

  if (grepl("hg19", type)) {

    path <- c(hg19_OMNI   = "interpolated_OMNI",
              hg19_hapmap = "interpolated_from_hapmap")[type]

    bigparallelr::register_parallel(ncores)

    foreach(chr = 1:22, .combine = "rbind") %dopar% {

      basename <- paste0("chr", chr, `if`(type == "hg19_OMNI", ".OMNI", ""),
                         ".interpolated_genetic_map")
      mapfile <- file.path(dir, basename)

      if (!file.exists(mapfile)) {
        url <- paste0("https://github.com/joepickrell/1000-genomes-genetic-maps/",
                      "raw/master/", path, "/", basename, ".gz")
        gzfile <- paste0(mapfile, ".gz")
        utils::download.file(url, destfile = gzfile, quiet = TRUE)
        R.utils::gunzip(gzfile)
      }

      cbind(chr = chr,
            bigreadr::fread2(mapfile, showProgress = FALSE, nThread = 1,
                             select = 2:3, col.names = c("pos", "pos_cM")))
    }

  } else if (type == "hg38_price") {

    mapfile <- file.path(dir, "genetic_map_hg38_withX.txt")

    if (!file.exists(mapfile)) {
      url <- "https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"
      gzfile <- paste0(mapfile, ".gz")
      utils::download.file(url, destfile = gzfile, quiet = TRUE)
      R.utils::gunzip(gzfile)
    }

    bigreadr::fread2(mapfile, showProgress = FALSE, nThread = ncores,
                     select = c(1, 2, 4), col.names = c("chr", "pos", "pos_cM"))

  }

}

################################################################################

#' Interpolate to genetic positions
#'
#' This function uses linear interpolation, whereas `snp_asGeneticPos()` uses
#' nearest neighbors.
#'
#' @inheritParams bigsnpr-package
#' @param genetic_map A data frame with 3 columns: `chr`, `pos`, and `pos_cM`.
#'   You can get it using [download_genetic_map()].
#'
#' @return The new vector of genetic positions.
#' @export
#'
snp_asGeneticPos2 <- function(infos.chr, infos.pos, genetic_map) {

  assert_lengths(infos.chr, infos.pos)
  assert_df_with_names(genetic_map, c("chr", "pos", "pos_cM"))

  genetic_map <- split(genetic_map, genetic_map$chr)

  new_pos <- rep(NA_real_, length(infos.chr))

  ind_chr <- split(seq_along(infos.chr), infos.chr)
  for (chr in names(ind_chr)) {
    ind.chr <- ind_chr[[chr]]
    ref <- genetic_map[[chr]]
    if (is.null(ref)) stop2("Chromosome '%s' not found in `genetic_map`.", chr)
    keep_unique <- !duplicated(ref$pos)
    new_pos[ind.chr] <- stats::approx(
      x = ref$pos[keep_unique], y = ref$pos_cM[keep_unique],
      xout = infos.pos[ind.chr], rule = 2)$y
  }

  new_pos
}

################################################################################
