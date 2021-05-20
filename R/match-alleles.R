################################################################################

flip_strand <- function(allele) {

  assert_package("dplyr")

  dplyr::case_when(
    allele == "A" ~ "T",
    allele == "C" ~ "G",
    allele == "T" ~ "A",
    allele == "G" ~ "C",
    TRUE ~ NA_character_
  )
}

#' Match alleles
#'
#' Match alleles between summary statistics and SNP information.
#' Match by ("chr", "a0", "a1") and ("pos" or "rsid"), accounting for possible
#' strand flips and reverse reference alleles (opposite effects).
#'
#' @param sumstats A data frame with columns "chr", "pos", "a0", "a1" and "beta".
#' @param info_snp A data frame with columns "chr", "pos", "a0" and "a1".
#' @param strand_flip Whether to try to flip strand? (default is `TRUE`)
#'   If so, ambiguous alleles A/T and C/G are removed.
#' @param join_by_pos Whether to join by chromosome and position (default),
#'   or instead by rsid.
#' @param remove_dups Whether to remove duplicates (same physical position)?
#'   Default is `TRUE`.
#' @param match.min.prop Minimum proportion of variants in the smallest data
#'   to be matched, otherwise stops with an error. Default is `50%`.
#'
#' @return A single data frame with matched variants. Values in column `$beta`
#'   are multiplied by -1 for variants with alleles reversed.
#' @export
#'
#' @seealso [snp_modifyBuild]
#'
#' @import data.table
#'
#' @example examples/example-match.R
snp_match <- function(sumstats, info_snp,
                      strand_flip = TRUE,
                      join_by_pos = TRUE,
                      remove_dups = TRUE,
                      match.min.prop = 0.5) {

  sumstats <- as.data.frame(sumstats)
  info_snp <- as.data.frame(info_snp)

  sumstats$`_NUM_ID_` <- rows_along(sumstats)
  info_snp$`_NUM_ID_` <- rows_along(info_snp)

  min_match <- match.min.prop * min(nrow(sumstats), nrow(info_snp))

  join_by <- c("chr", NA, "a0", "a1")
  join_by[2] <- `if`(join_by_pos, "pos", "rsid")

  if (!all(c(join_by, "beta") %in% names(sumstats)))
    stop2("Please use proper names for variables in 'sumstats'. Expected '%s'.",
          paste(c(join_by, "beta"), collapse = ", "))
  if (!all(c(join_by, "pos") %in% names(info_snp)))
    stop2("Please use proper names for variables in 'info_snp'. Expected '%s'.",
          paste(unique(c(join_by, "pos")), collapse = ", "))

  message2("%s variants to be matched.", format(nrow(sumstats), big.mark = ","))

  # first filter to fasten
  sumstats <- sumstats[vctrs::vec_in(sumstats[, join_by[1:2]],
                                     info_snp[, join_by[1:2]]), ]
  if (nrow(sumstats) == 0)
    stop2("No variant has been matched.")

  # augment dataset to match reverse alleles
  if (strand_flip) {
    is_ambiguous <- with(sumstats, paste(a0, a1) %in% c("A T", "T A", "C G", "G C"))
    message2("%s ambiguous SNPs have been removed.",
             format(sum(is_ambiguous), big.mark = ","))
    sumstats2 <- sumstats[!is_ambiguous, ]
    sumstats3 <- sumstats2
    sumstats2$`_FLIP_` <- FALSE
    sumstats3$`_FLIP_` <- TRUE
    sumstats3$a0 <- flip_strand(sumstats2$a0)
    sumstats3$a1 <- flip_strand(sumstats2$a1)
    sumstats3 <- rbind(sumstats2, sumstats3)
  } else {
    sumstats3 <- sumstats
    sumstats3$`_FLIP_` <- FALSE
  }

  sumstats4 <- sumstats3
  sumstats3$`_REV_` <- FALSE
  sumstats4$`_REV_` <- TRUE
  sumstats4$a0 <- sumstats3$a1
  sumstats4$a1 <- sumstats3$a0
  sumstats4$beta <- -sumstats3$beta
  sumstats4 <- rbind(sumstats3, sumstats4)

  matched <- merge(as.data.table(sumstats4), as.data.table(info_snp),
                   by = join_by, all = FALSE, suffixes = c(".ss", ""))

  if (remove_dups) {
    dups <- vctrs::vec_duplicate_detect(matched[, c("chr", "pos")])
    if (any(dups)) {
      matched <- matched[!dups, ]
      message2("Some duplicates were removed.")
    }
  }

  message2("%s variants have been matched; %s were flipped and %s were reversed.",
           format(nrow(matched),         big.mark = ","),
           format(sum(matched$`_FLIP_`), big.mark = ","),
           format(sum(matched$`_REV_`),  big.mark = ","))

  if (nrow(matched) < min_match)
    stop2("Not enough variants have been matched.")

  as.data.frame(matched[, c("_FLIP_", "_REV_") := NULL][order(chr, pos)])
}

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
#'
#' @references
#' Hinrichs, Angela S., et al. "The UCSC genome browser database: update 2006."
#' Nucleic acids research 34.suppl_1 (2006): D590-D598.
#'
#' @return Input data frame `info_snp` with column "pos" in the new build.
#' @export
#'
snp_modifyBuild <- function(info_snp, liftOver, from = "hg18", to = "hg19") {

  if (!all(c("chr", "pos") %in% names(info_snp)))
    stop2("Please use proper names for variables in 'info_snp'. Expected %s.",
          "'chr' and 'pos'")

  # Make sure liftOver is executable
  liftOver <- normalizePath(liftOver)
  make_executable(liftOver)

  # Need BED UCSC file for liftOver
  BED <- tempfile(fileext = ".BED")
  info_BED <- with(info_snp, data.frame(
    paste0("chr", chr), pos0 = pos - 1L, pos, id = rows_along(info_snp)))
  bigreadr::fwrite2(info_BED, BED, col.names = FALSE, sep = " ")

  # Need chain file
  url <- paste0("ftp://hgdownload.cse.ucsc.edu/goldenPath/", from, "/liftOver/",
                from, "To", tools::toTitleCase(to), ".over.chain.gz")
  chain <- tempfile(fileext = ".over.chain.gz")
  utils::download.file(url, destfile = chain)

  # Run liftOver (usage: liftOver oldFile map.chain newFile unMapped)
  lifted <- tempfile(fileext = ".BED")
  unmapped <- tempfile(fileext = ".txt")
  system(paste(liftOver, BED, chain, lifted, unmapped))

  # readLines(lifted, n = 5)
  new_pos <- bigreadr::fread2(lifted)

  # readLines(unmapped, n = 6)
  bad <- grep("^#", readLines(unmapped), value = TRUE, invert = TRUE)
  message2("%d variants have not been mapped.", length(bad))

  info_snp$pos <- NA
  info_snp$pos[new_pos$V4] <- new_pos$V3
  info_snp
}

################################################################################

#' Determine reference divergence
#'
#' Determine reference divergence while accounting for strand flips.
#' **This does not remove ambiguous alleles.**
#'
#' @param ref1 The reference alleles of the first dataset.
#' @param alt1 The alternative alleles of the first dataset.
#' @param ref2 The reference alleles of the second dataset.
#' @param alt2 The alternative alleles of the second dataset.
#'
#' @return A logical vector whether the references alleles are the same.
#'   Missing values can result from missing values in the inputs or from
#'   ambiguous matching (e.g. matching A/C and A/G).
#' @export
#'
#' @seealso [snp_match()]
#'
#' @examples
#' same_ref(ref1 = c("A", "C", "T", "G", NA),
#'          alt1 = c("C", "T", "C", "A", "A"),
#'          ref2 = c("A", "C", "A", "A", "C"),
#'          alt2 = c("C", "G", "G", "G", "A"))
same_ref <- function(ref1, alt1, ref2, alt2) {

  # ACTG <- c("A", "C", "T", "G")
  # REV_ACTG <- stats::setNames(c("T", "G", "A", "C"), ACTG)
  #
  # decoder <- expand.grid(list(ACTG, ACTG, ACTG, ACTG)) %>%
  #   dplyr::mutate(status = dplyr::case_when(
  #     # BAD: same reference/alternative alleles in a dataset
  #     (Var1 == Var2) | (Var3 == Var4) ~ NA,
  #     # GOOD/TRUE: same reference/alternative alleles between datasets
  #     (Var1 == Var3) & (Var2 == Var4) ~ TRUE,
  #     # GOOD/FALSE: reverse reference/alternative alleles
  #     (Var1 == Var4) & (Var2 == Var3) ~ FALSE,
  #     # GOOD/TRUE: same reference/alternative alleles after strand flip
  #     (REV_ACTG[Var1] == Var3) & (REV_ACTG[Var2] == Var4) ~ TRUE,
  #     # GOOD/FALSE: reverse reference/alternative alleles after strand flip
  #     (REV_ACTG[Var1] == Var4) & (REV_ACTG[Var2] == Var3) ~ FALSE,
  #     # BAD: the rest
  #     TRUE ~ NA
  #   )) %>%
  #   reshape2::acast(Var1 ~ Var2 ~ Var3 ~ Var4, value.var = "status")

  decoder <- structure(
    c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE,
      NA, NA, FALSE, NA, NA, NA, NA, NA, NA, TRUE, NA, NA, FALSE, NA, NA, NA,
      TRUE, NA, NA, NA, NA, NA, FALSE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
      TRUE, NA, NA, FALSE, NA, NA, TRUE, NA, NA, FALSE, NA, NA, NA, NA, FALSE,
      NA, NA, TRUE, NA, NA, NA, NA, NA, NA, FALSE, NA, NA, TRUE, NA, NA, NA, NA,
      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, FALSE, NA,
      NA, TRUE, NA, NA, FALSE, NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA,
      NA, TRUE, NA, NA, NA, NA, NA, FALSE, NA, NA, NA, NA, FALSE, NA, NA, NA,
      NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, NA, NA, FALSE,
      NA, NA, TRUE, NA, NA, FALSE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
      NA, NA, NA, NA, NA, NA, NA, NA, NA, TRUE, NA, NA, FALSE, NA, NA, NA, NA,
      NA, NA, TRUE, NA, NA, FALSE, NA, NA, NA, NA, FALSE, NA, NA, TRUE, NA, NA,
      FALSE, NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, FALSE, NA,
      NA, NA, NA, NA, TRUE, NA, NA, NA, FALSE, NA, NA, TRUE, NA, NA, NA, NA, NA,
      NA, FALSE, NA, NA, TRUE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
      NA, NA, NA, NA, NA),
    .Dim = rep(4, 4), .Dimnames = rep(list(c("A", "C", "T", "G")), 4)
  )

  to_decode <- do.call("cbind", lapply(list(ref1, alt1, ref2, alt2), as.character))
  decoder[to_decode]
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
#'   Default is `tempdir()`. Directly use files there if already present.
#' @param rsid If providing rsIDs, the matching is performed using those
#'   (instead of positions) and variants not matched are interpolated using
#'   spline interpolation of variants that have been matched.
#'
#' @return The new vector of genetic positions.
#' @export
#'
snp_asGeneticPos <- function(infos.chr, infos.pos, dir = tempdir(), ncores = 1,
                             rsid = NULL) {

  assert_package("R.utils")
  assert_lengths(infos.chr, infos.pos)
  if (!is.null(rsid)) assert_lengths(rsid, infos.pos)

  snp_split(infos.chr, function(ind.chr, pos, dir, rsid) {

    chr <- attr(ind.chr, "chr")
    basename <- paste0("chr", chr, ".OMNI.interpolated_genetic_map")
    mapfile <- file.path(dir, basename)
    if (!file.exists(mapfile)) {
      url <- paste0("https://github.com/joepickrell/1000-genomes-genetic-maps/",
                    "raw/master/interpolated_OMNI/", basename, ".gz")
      gzfile <- paste0(mapfile, ".gz")
      utils::download.file(url, destfile = gzfile, quiet = TRUE)
      R.utils::gunzip(gzfile)
    }
    map.chr <- bigreadr::fread2(mapfile, showProgress = FALSE)

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
          stats::spline(pos.chr, new_pos, xout = pos.chr[indNA], method = "hyman", ties = "ordered")$y)
      }
    }

    new_pos

  }, combine = "c", pos = infos.pos, dir = dir, rsid = rsid, ncores = ncores)
}

################################################################################
