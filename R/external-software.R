################################################################################

# https://github.com/r-lib/rappdirs/blob/master/R/utils.r
get_os <- function() {
  if (.Platform$OS.type == "windows") {
    "Windows"
  } else if (Sys.info()[["sysname"]] == "Darwin") {
    "Mac"
  } else if (.Platform$OS.type == "unix") {
    "Unix"
  } else {
    stop("Unknown OS")
  }
}

#' Download PLINK 1.9
#'
#' Download PLINK 1.9 from \url{https://www.cog-genomics.org/plink2}.
#'
#' @param dir The directory where to put the PLINK executable.
#'   Default is a temporary directory.
#'
#' @return The path of the downloaded PLINK executable.
#'
#' @export
#'
download_plink <- function(dir = tempdir()) {

  myOS <- get_os()
  PLINK <- file.path(dir, `if`(myOS == "Windows", "plink.exe", "plink"))
  if (file.exists(PLINK)) return(PLINK)

  # https://regex101.com/r/jC8nB0/143
  plink.names  <- gsubfn::strapply(
    X = readLines("https://www.cog-genomics.org/plink2"),
    pattern = "(/static/bin/.*/plink_.+?(?<!dev)\\.zip)",
    simplify = "c",
    perl = TRUE
  )
  plink.builds <- data.frame(
    url = paste0("https://www.cog-genomics.org", plink.names),
    OS = c(rep("Unix", 2), "Mac", rep("Windows", 2)),
    arch = c(64, 32, 64, 64, 32),
    stringsAsFactors = FALSE
  )

  myArch <- 8 * .Machine$sizeof.pointer
  url <- subset(plink.builds, OS == myOS & arch == myArch)[["url"]]

  utils::download.file(url, destfile = (plink.zip <- tempfile(fileext = ".zip")))
  PLINK <- utils::unzip(plink.zip,
                        files = basename(PLINK),
                        exdir = dirname(PLINK))
  Sys.chmod(PLINK, mode = (file.info(PLINK)$mode | "111"))

  PLINK
}

################################################################################

#' Download Beagle 4.1
#'
#' Download Beagle 4.1 from
#' \url{https://faculty.washington.edu/browning/beagle/beagle.html}
#'
#' @param dir The directory where to put the Beagle Java Archive.
#'   Default is a temporary directory.
#'
#' @return The path of the downloaded Beagle Java Archive.
#'
#' @export
#'
download_beagle <- function(dir = tempdir()) {

  url <- "https://faculty.washington.edu/browning/beagle/"

  # https://regex101.com/r/jC8nB0/141
  jar  <- gsubfn::strapply(
    X = readLines(paste0(url, "beagle.html")),
    pattern = "(beagle.+?\\.jar)",
    simplify = "c",
    perl = TRUE
  )[[1]]

  dest <- file.path(dir, jar)

  if (!file.exists(dest)) {
    utils::download.file(paste0(url, jar), destfile = dest)
    Sys.chmod(dest, mode = (file.info(dest)$mode | "111"))
  }

  dest
}

################################################################################

#' Quality Control
#'
#' Quality Control (QC) and possible conversion to *bed*/*bim*/*fam* files
#' using [**PLINK 1.9**](https://www.cog-genomics.org/plink2).
#'
#' @param plink.path Path to the executable of PLINK 1.9.
#' @param prefix.in Prefix (path without extension) of the dataset to be QCed.
#' @param file.type Type of the dataset to be QCed. Default is `"--bfile"` and
#'   corresponds to bed/bim/fam files. You can also use `"--file"` for ped/map
#'   files or `"--vcf"` for a VCF file. More information can be found at
#'   \url{https://www.cog-genomics.org/plink/1.9/input}.
#' @param prefix.out Prefix (path without extension) of the bed/bim/fam dataset
#'   to be created. Default is created by appending `"_QC"` to `prefix.in`.
#' @param maf Minimum Minor Allele Frequency (MAF) for a SNP to be kept.
#'   Default is `0.01`.
#' @param geno Maximum proportion of missing values for a SNP to be kept.
#'   Default is `0.1`.
#' @param mind Maximum proportion of missing values for a sample to be kept.
#'   Default is `0.1`.
#' @param hwe Filters out all variants which have Hardy-Weinberg equilibrium
#'   exact test p-value below the provided threshold. Default is `1e-50`.
#' @param autosome.only Whether to exclude all unplaced and non-autosomal
#'   variants? Default is `FALSE`.
#' @param extra.options Other options to be passed to PLINK as a string. More
#'   options can be found at \url{https://www.cog-genomics.org/plink2/filter}.
#'
#' @return The path of the newly created bedfile.
#'
#' @export
#'
#' @references
#' Purcell, Shaun, Benjamin Neale, Kathe Todd-Brown, Lori Thomas,
#' Manuel A R Ferreira, David Bender, Julian Maller, et al. 2007.
#' *PLINK: a tool set for whole-genome association and population-based linkage
#' analyses.* American Journal of Human Genetics 81 (3). Elsevier: 559–75.
#' \url{http://dx.doi.org/10.1086/519795}.
#'
#' Chang, Christopher C, Carson C Chow, Laurent CAM Tellier,
#' Shashaank Vattikuti, Shaun M Purcell, and James J Lee. 2015.
#' *Second-generation PLINK: rising to the challenge of larger and richer
#' datasets.* GigaScience 4 (1): 7.
#' \url{http://dx.doi.org/10.1186/s13742-015-0047-8}.
#'
#' @seealso [download_plink] [snp_plinkIBDQC]
#'
#' @examples
#' bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#' prefix  <- sub("\\.bed$", "", bedfile)
#' plink <- download_plink()
#' test <- snp_plinkQC(plink.path = plink,
#'                     prefix.in = prefix,
#'                     prefix.out = tempfile(),
#'                     file.type = "--bfile",  # the default (for ".bed")
#'                     maf = 0.05,
#'                     geno = 0.05,
#'                     mind = 0.05,
#'                     hwe = 1e-10,
#'                     autosome.only = TRUE)
#' test
#'
snp_plinkQC <- function(plink.path,
                        prefix.in,
                        file.type = "--bfile",
                        prefix.out = paste0(prefix.in, "_QC"),
                        maf = 0.01,
                        geno = 0.1,
                        mind = 0.1,
                        hwe = 1e-50,
                        autosome.only = FALSE,
                        extra.options = "") {

  # new bedfile, check if already exists
  bedfile.out <- paste0(prefix.out, ".bed")
  assert_noexist(bedfile.out)

  # --vcf expects a complete filename
  if (file.type == "--vcf") prefix.in <- paste0(prefix.in, ".vcf")

  # call PLINK 1.9
  system(
    paste(
      plink.path,
      file.type, prefix.in,
      "--maf", maf,
      "--mind", mind,
      "--geno", geno,
      "--hwe", hwe,
      `if`(autosome.only, "--autosome", ""),
      "--make-bed",
      "--out", prefix.out,
      extra.options
    )
  )

  # return path to new bedfile
  bedfile.out
}

################################################################################

#' Remove samples
#'
#' Create new *bed*/*bim*/*fam* files by removing samples with PLINK.
#'
#' @inheritParams snp_plinkQC
#' @param bedfile.in Path to the input bedfile.
#' @param bedfile.out Path to the output bedfile.
#' @inheritParams snp_getSampleInfos
#'
#' @seealso [download_plink]
#'
#' @return The path of the new bedfile.
#' @export
#'
snp_plinkRmSamples <- function(plink.path,
                               bedfile.in,
                               bedfile.out,
                               df.or.files,
                               col.family.ID = 1,
                               col.sample.ID = 2,
                               ...) {

  assert_noexist(bedfile.out)

  if (is.data.frame(df.or.files)) {
    data.infos <- df.or.files
  } else if (is.character(df.or.files)) {
    data.infos <- foreach(f = df.or.files, .combine = 'rbind') %do% {
      data.table::fread(f, data.table = FALSE, ...)
    }
  } else {
    stop2("'df.or.files' must be a data.frame or a vector of file paths.")
  }

  tmpfile <- tempfile()
  write.table2(data.infos[, c(col.family.ID, col.sample.ID)], tmpfile)
  # Make new bed with extraction
  system(
    paste(
      plink.path,
      "--bfile", sub("\\.bed$", "", bedfile.in),
      "--make-bed",
      "--out", sub("\\.bed$", "", bedfile.out),
      "--remove", tmpfile
    )
  )

  bedfile.out
}

################################################################################

#' Identity-by-descent
#'
#' Quality Control based on Identity-by-descent (IBD) computed by
#' [**PLINK 1.9**](https://www.cog-genomics.org/plink2)
#' using its method-of-moments.
#'
#' @inheritParams snp_plinkRmSamples
#' @param bedfile.out Path to the output bedfile. Default is created by
#'   appending `"_norel"` to `prefix.in` (`bedfile.in` without extension).
#' @param pi.hat PI_HAT value threshold for individuals (first by pairs)
#'   to be excluded. Default is `0.08`.
#' @inheritParams bigsnpr-package
#' @param pruning.args A vector of 2 pruning parameters, respectively
#'   the window size (in variant count) and the pairwise $r^2$ threshold
#'   (the step size is fixed to 1). Default is `c(100, 0.2)`.
#' @param extra.options Other options to be passed to PLINK as a string
#'   (for the IBD part). More options can be found at
#'   \url{https://www.cog-genomics.org/plink/1.9/ibd}.
#' @param do.blind.QC Whether to do QC with `pi.hat` without visual inspection.
#'   Default is `TRUE`. If `FALSE`, return the `data.frame` of the corresponding
#'   ".genome" file without doing QC. One could use
#'   `ggplot2::qplot(Z0, Z1, data = mydf, col = RT)` for visual inspection.
#'
#' @return The path of the new bedfile.
#'   If no sample is filter, no new bed/bim/fam files are created and
#'   then the path of the input bedfile is returned.
#' @export
#'
#' @inherit snp_plinkQC references
#'
#' @seealso [download_plink] [snp_plinkQC]
#'
#' @examples
#' bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#' plink <- download_plink()
#' test <- snp_plinkIBDQC(plink, bedfile,
#'                        bedfile.out = tempfile(fileext = ".bed"),
#'                        ncores = 2)
#' test
#' test2 <- snp_plinkIBDQC(plink, bedfile,
#'                         do.blind.QC = FALSE,
#'                         ncores = 2)
#' str(test2)
#' library(ggplot2)
#' qplot(Z0, Z1, data = test2, col = RT)
#' qplot(y = PI_HAT, data = test2) +
#'   geom_hline(yintercept = 0.2, color = "blue", linetype = 2)
#' snp_plinkRmSamples(plink, bedfile,
#'                    bedfile.out = tempfile(fileext = ".bed"),
#'                    df.or.files = subset(test2, PI_HAT > 0.2))
#'
snp_plinkIBDQC <- function(plink.path,
                           bedfile.in,
                           bedfile.out = NULL,
                           pi.hat = 0.08,
                           ncores = 1,
                           pruning.args = c(100, 0.2),
                           do.blind.QC = TRUE,
                           extra.options = "") {

  # temporary file
  tmpfile <- tempfile()

  # check extension of file
  assert_ext(bedfile.in, "bed")
  # get file without extension
  prefix.in <- sub("\\.bed$", "", bedfile.in)

  # get possibly new file
  if (is.null(bedfile.out)) bedfile.out <- paste0(prefix.in, "_norel.bed")
  if (do.blind.QC) assert_noexist(bedfile.out)

  # prune if desired
  if (!is.null(pruning.args)) {
    system(
      paste(
        plink.path,
        "--bfile", prefix.in,
        "--indep-pairwise", pruning.args[1], 1, pruning.args[2],
        "--out", tmpfile
      )
    )
    opt.pruning <- paste0("--extract ", tmpfile, ".prune.in")
  } else {
    opt.pruning <- ""
  }

  # compute IBD
  system(
    paste(
      plink.path,
      "--bfile", prefix.in,
      opt.pruning,
      "--genome",
      "--min", pi.hat,
      "--out", tmpfile,
      "--threads", ncores,
      extra.options
    )
  )

  # get genomefile as a data.frame
  tmp <- data.table::fread(paste0(tmpfile, ".genome"), data.table = FALSE)
  if (nrow(tmp)) { # if there are samples to filter

    if (do.blind.QC) {
      snp_plinkRmSamples(plink.path, bedfile.in, bedfile.out, tmp)
    } else {
      tmp
    }

  } else {

    message2("No pair of samples has PI_HAT > %s.", pi.hat)
    message("Returning input bedfile..")
    bedfile.in

  }
}

################################################################################

#' Imputation
#'
#' Imputation using **Beagle** version 4.
#'
#' Downloads and more information can be found at the following websites
#' - [PLINK](https://www.cog-genomics.org/plink2),
#' - [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html).
#'
#' @param beagle.path Path to the executable of Beagle v4+.
#' @inheritParams snp_plinkRmSamples
#' @param bedfile.out Path to the output bedfile. Default is created by
#'   appending `"_impute"` to `prefix.in` (`bedfile.in` without extension).
#' @param memory.max Max memory (in GB) to be used. It is internally rounded
#'   to be an integer. Default is `3`.
#' @param ncores Number of cores to be used. Default is `1`. An usually good
#'   value for this parameter is `ncores = parallel::detectCores() - 1`.
#' @param extra.options Other options to be passed to Beagle as a string. More
#'   options can be found at Beagle's website.
#'
#' @references B L Browning and S R Browning (2016).
#' Genotype imputation with millions of reference samples.
#' Am J Hum Genet 98:116-126.
#' \url{dx.doi.org/doi:10.1016/j.ajhg.2015.11.020}
#'
#' @seealso [download_plink]
#'
#' @return The path of the new bedfile.
#' @export
snp_beagleImpute <- function(beagle.path,
                             plink.path,
                             bedfile.in,
                             bedfile.out = NULL,
                             memory.max = 3,
                             ncores = 1,
                             extra.options = "") {

  # get input and output right
  prefix.in <- sub("\\.bed$", "", bedfile.in)
  if (is.null(bedfile.out)) bedfile.out <- paste0(prefix.in, "_impute.bed")
  assert_noexist(bedfile.out)
  prefix.out <- sub("\\.bed$", "", bedfile.out)

  # get temporary files
  tmpfile1 <- tempfile()
  tmpfile2 <- tempfile()

  # Convert bed/bim/fam to vcf
  system(
    paste(
      plink.path,
      "--bfile", prefix.in,
      "--recode vcf bgz", # .vcf.gz extension
      "--out", tmpfile1,
      "--threads", ncores
    )
  )
  vcf1 <- paste0(tmpfile1, ".vcf.gz")
  on.exit(file.remove(vcf1), add = TRUE)

  # Impute vcf with Beagle version 4
  system(
    paste(
      "java", sprintf("-Xmx%dg", round(memory.max)),
      "-jar", beagle.path,
      paste0("gt=", vcf1),
      paste0("out=", tmpfile2),
      paste0("nthreads=", ncores),
      extra.options
    )
  )
  vcf2 <- paste0(tmpfile2, ".vcf.gz")
  on.exit(file.remove(vcf2), add = TRUE)

  # Convert back vcf to bed/bim/fam
  system(
    paste(
      plink.path,
      "--vcf", vcf2,
      "--out", prefix.out
    )
  )

  # return path to new bedfile
  bedfile.out
}

################################################################################
