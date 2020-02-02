################################################################################

make_executable <- function(exe) {
  Sys.chmod(exe, mode = (file.info(exe)$mode | "111"))
}

system_verbose <- function(..., verbose) {
  system(..., ignore.stdout = !verbose, ignore.stderr = !verbose)
}

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

################################################################################

# Thanks @richfitz for this
get_pattern <- function(x, pattern) {
  sub(
    pattern = pattern,
    replacement = "\\1",
    x = grep(pattern, x, value = TRUE)
  )
}

################################################################################

#' Download PLINK
#'
#' Download PLINK 1.9 from \url{http://www.cog-genomics.org/plink2}.
#'
#' @param dir The directory where to put the PLINK executable.
#'   Default is a temporary directory.
#' @param overwrite Whether to overwrite file? Default is `FALSE`.
#' @param verbose Whether to output details of downloading. Default is `TRUE`.
#'
#' @return The path of the downloaded PLINK executable.
#'
#' @export
#'
download_plink <- function(dir = tempdir(), overwrite = FALSE, verbose = TRUE) {

  myOS <- get_os()
  PLINK <- file.path(dir, `if`(myOS == "Windows", "plink.exe", "plink"))
  if (!overwrite && file.exists(PLINK)) return(PLINK)

  plink.names <- get_pattern(
    x = readLines("http://www.cog-genomics.org/plink2"),
    # http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190304.zip
    pattern = ".*(http://s3.amazonaws.com/plink1-assets/plink_.+?\\.zip).*"
  )
  plink.builds <- data.frame(
    url = plink.names,
    OS = c(rep("Unix", 2), "Mac", rep("Windows", 2)),
    arch = c(64, 32, 64, 64, 32),
    stringsAsFactors = FALSE
  )

  myArch <- 8 * .Machine$sizeof.pointer
  url <- subset(plink.builds, OS == myOS & arch == myArch)[["url"]]

  plink.zip <- tempfile(fileext = ".zip")
  utils::download.file(url, destfile = plink.zip, quiet = !verbose)
  PLINK <- utils::unzip(plink.zip,
                        files = basename(PLINK),
                        exdir = dirname(PLINK))
  make_executable(PLINK)

  PLINK
}

################################################################################

#' Download PLINK
#'
#' Download PLINK 2.0 from \url{http://www.cog-genomics.org/plink/2.0/}.
#'
#' @param AVX2 Whether to download the AVX2 version? This is only available for
#'   64 bits architectures. Default tries to detect if system supports AVX2.
#'
#' @export
#'
#' @rdname download_plink
#'
download_plink2 <- function(dir = tempdir(), AVX2 = has_avx2(),
                            overwrite = FALSE, verbose = TRUE) {

  myOS <- get_os()
  PLINK2 <- file.path(dir, `if`(myOS == "Windows", "plink2.exe", "plink2"))
  if (!overwrite && file.exists(PLINK2)) return(PLINK2)

  plink.names <- get_pattern(
    x = readLines("http://www.cog-genomics.org/plink/2.0/"),
    # http://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20190527.zip
    pattern = ".*(http://s3.amazonaws.com/plink2-assets/plink2_.+?\\.zip).*"
  )
  plink.builds <- data.frame(
    url  = plink.names,
    OS   = c(rep("Unix", 3),      rep("Mac", 2),  rep("Windows", 3)),
    arch = c(64, 64, 32,          64, 64,         64, 64, 32),
    avx2 = c(TRUE, FALSE, FALSE,  TRUE, FALSE,    TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  myArch <- 8 * .Machine$sizeof.pointer
  if (myArch == 32) AVX2 <- FALSE
  url <- subset(plink.builds, OS == myOS & arch == myArch & avx2 == AVX2)[["url"]]

  plink.zip <- tempfile(fileext = ".zip")
  utils::download.file(url, destfile = plink.zip, quiet = !verbose)
  PLINK2 <- utils::unzip(plink.zip,
                         files = basename(PLINK2),
                         exdir = dirname(PLINK2))
  make_executable(PLINK2)

  PLINK2
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

  jar <- get_pattern(
    x = readLines(paste0(url, "beagle.html")),
    pattern = ".*(beagle.+?\\.jar).*"
  )

  beagle <- file.path(dir, jar)
  if (!file.exists(beagle)) {
    utils::download.file(paste0(url, jar), destfile = beagle)
    make_executable(beagle)
  }

  beagle
}

################################################################################

#' Quality Control
#'
#' Quality Control (QC) and possible conversion to *bed*/*bim*/*fam* files
#' using [**PLINK 1.9**](http://www.cog-genomics.org/plink2).
#'
#' @param plink.path Path to the executable of PLINK 1.9.
#' @param prefix.in Prefix (path without extension) of the dataset to be QCed.
#' @param file.type Type of the dataset to be QCed. Default is `"--bfile"` and
#'   corresponds to bed/bim/fam files. You can also use `"--file"` for ped/map
#'   files or `"--vcf"` for a VCF file. More information can be found at
#'   \url{http://www.cog-genomics.org/plink/1.9/input}.
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
#'   options can be found at \url{http://www.cog-genomics.org/plink2/filter}.
#'   If using PLINK 2.0, you could e.g. use `"--king-cutoff 0.0884"` to remove
#'   some related samples at the same time of quality controls.
#' @param verbose Whether to show PLINK log? Default is `TRUE`.
#'
#' @return The path of the newly created bedfile.
#'
#' @export
#'
#' @references
#' Chang, Christopher C, Carson C Chow, Laurent CAM Tellier,
#' Shashaank Vattikuti, Shaun M Purcell, and James J Lee. 2015.
#' *Second-generation PLINK: rising to the challenge of larger and richer
#' datasets.* GigaScience 4 (1): 7.
#' \url{http://dx.doi.org/10.1186/s13742-015-0047-8}.
#'
#' @seealso [download_plink] [snp_plinkIBDQC]
#'
#' @examples
#' \dontrun{
#'
#' bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#' prefix  <- sub_bed(bedfile)
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
#' }
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
                        extra.options = "",
                        verbose = TRUE) {

  # new bedfile, check if already exists
  bedfile.out <- paste0(prefix.out, ".bed")
  assert_noexist(bedfile.out)

  # --vcf expects a complete filename
  if (file.type == "--vcf") prefix.in <- paste0(prefix.in, ".vcf")

  # call PLINK 1.9
  system_verbose(
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
    ),
    verbose = verbose
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
                               ...,
                               verbose = TRUE) {

  assert_noexist(bedfile.out)

  if (is.data.frame(df.or.files)) {
    data.infos <- df.or.files
  } else if (is.character(df.or.files)) {
    data.infos <- bigreadr::fread2(df.or.files, ...)
  } else {
    stop2("'df.or.files' must be a data.frame or a vector of file paths.")
  }

  tmpfile <- tempfile()
  write.table2(data.infos[, c(col.family.ID, col.sample.ID)], tmpfile)
  # Make new bed with extraction
  system_verbose(
    paste(
      plink.path,
      "--bfile", sub_bed(bedfile.in),
      "--make-bed",
      "--out", sub_bed(bedfile.out),
      "--remove", tmpfile
    ),
    verbose = verbose
  )

  bedfile.out
}

################################################################################

#' Identity-by-descent
#'
#' Quality Control based on Identity-by-descent (IBD) computed by
#' [**PLINK 1.9**](http://www.cog-genomics.org/plink2)
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
#'   \url{http://www.cog-genomics.org/plink/1.9/ibd}.
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
#' @seealso [download_plink] [snp_plinkQC] [snp_plinkKINGQC]
#'
#' @examples
#' \dontrun{
#'
#' bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#' plink <- download_plink()
#'
#' bedfile <- snp_plinkIBDQC(plink, bedfile,
#'                           bedfile.out = tempfile(fileext = ".bed"),
#'                           ncores = 2)
#'
#' df_rel <- snp_plinkIBDQC(plink, bedfile, do.blind.QC = FALSE, ncores = 2)
#' str(df_rel)
#'
#' library(ggplot2)
#' qplot(Z0, Z1, data = df_rel, col = RT)
#' qplot(y = PI_HAT, data = df_rel) +
#'   geom_hline(yintercept = 0.2, color = "blue", linetype = 2)
#' snp_plinkRmSamples(plink, bedfile,
#'                    bedfile.out = tempfile(fileext = ".bed"),
#'                    df.or.files = subset(df_rel, PI_HAT > 0.2))
#' }
#'
snp_plinkIBDQC <- function(plink.path,
                           bedfile.in,
                           bedfile.out = NULL,
                           pi.hat = 0.08,
                           ncores = 1,
                           pruning.args = c(100, 0.2),
                           do.blind.QC = TRUE,
                           extra.options = "",
                           verbose = TRUE) {

  # temporary file
  tmpfile <- tempfile()

  # get file without extension
  prefix.in <- sub_bed(bedfile.in)

  # get possibly new file
  if (is.null(bedfile.out)) bedfile.out <- paste0(prefix.in, "_norel.bed")
  if (do.blind.QC) assert_noexist(bedfile.out)

  # prune if desired
  if (!is.null(pruning.args)) {
    system_verbose(
      paste(
        plink.path,
        "--bfile", prefix.in,
        "--indep-pairwise", pruning.args[1], 1, pruning.args[2],
        "--out", tmpfile
      ),
      verbose = verbose
    )
    opt.pruning <- paste0("--extract ", tmpfile, ".prune.in")
  } else {
    opt.pruning <- ""
  }

  # compute IBD
  system_verbose(
    paste(
      plink.path,
      "--bfile", prefix.in,
      opt.pruning,
      "--genome",
      "--min", pi.hat,
      "--out", tmpfile,
      "--threads", ncores,
      extra.options
    ),
    verbose = verbose
  )

  # get genomefile as a data.frame
  tmp <- bigreadr::fread2(paste0(tmpfile, ".genome"))
  if (nrow(tmp) > 0) { # if there are samples to filter

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

#' Relationship-based pruning
#'
#' Quality Control based on KING-robust kinship estimator. More information can
#' be found at \url{http://www.cog-genomics.org/plink/2.0/distance#king_cutoff}.
#'
#' @param plink2.path Path to the executable of PLINK 2.
#' @inheritParams snp_plinkIBDQC
#' @param thr.king  Note that KING kinship coefficients are scaled such that
#'   duplicate samples have kinship 0.5, not 1. First-degree relations
#'   (parent-child, full siblings) correspond to ~0.25, second-degree relations
#'   correspond to ~0.125, etc. It is conventional to use a cutoff of ~0.354
#'   (2^-1.5, the geometric mean of 0.5 and 0.25) to screen for monozygotic
#'   twins and duplicate samples, ~0.177 (2^-2.5) to remove first-degree
#'   relations as well, and ~0.0884 (2^-3.5, **default**) to remove
#'   second-degree relations as well, etc.
#' @param extra.options Other options to be passed to PLINK2 as a string.
#' @param make.bed Whether to create new bed/bim/fam files (default).
#'   Otherwise, returns a table with coefficients of related pairs.
#'
#' @return See parameter `make-bed`.
#' @export
#'
#' @inherit snp_plinkQC references
#' @references
#' Manichaikul, Ani, Josyf C. Mychaleckyj, Stephen S. Rich, Kathy Daly,
#' Michele Sale, and Wei-Min Chen. "Robust relationship inference in genome-wide
#' association studies." Bioinformatics 26, no. 22 (2010): 2867-2873.
#'
#' @seealso [download_plink2] [snp_plinkQC]
#'
#' @examples
#' \dontrun{
#'
#' bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#' plink2 <- download_plink2(AVX2 = FALSE)
#'
#' bedfile2 <- snp_plinkKINGQC(plink2, bedfile,
#'                             bedfile.out = tempfile(fileext = ".bed"),
#'                             ncores = 2)
#'
#' df_rel <- snp_plinkKINGQC(plink2, bedfile, make.bed = FALSE, ncores = 2)
#' str(df_rel)
#' }
#'
snp_plinkKINGQC <- function(plink2.path,
                            bedfile.in,
                            bedfile.out = NULL,
                            thr.king = 2^-3.5,
                            make.bed = TRUE,
                            ncores = 1,
                            extra.options = "",
                            verbose = TRUE) {

  # check PLINK version
  v <- system(paste(plink2.path, "--version"), intern = TRUE)
  if (substr(v, 1, 8) != "PLINK v2")
    stop2("This requires PLINK v2; got '%s' instead.", v)

  # get file without extension
  prefix.in <- sub_bed(bedfile.in)

  # get possibly new file
  if (make.bed) {

    if (is.null(bedfile.out)) bedfile.out <- paste0(prefix.in, "_norel.bed")
    assert_noexist(bedfile.out)

    # compute KING-robust kinship coefficients and filter
    system_verbose(
      paste(
        plink2.path,
        "--bfile", prefix.in,
        "--make-bed --king-cutoff", thr.king,
        "--out", sub_bed(bedfile.out),
        "--threads", ncores,
        extra.options
      ),
      verbose = verbose
    )

    bedfile.out

  } else {

    prefix.out <- tempfile()

    # compute table of KING-robust kinship coefficients
    system_verbose(
      paste(
        plink2.path,
        "--bfile", prefix.in,
        "--make-king-table --king-table-filter", thr.king,
        "--out", prefix.out,
        "--threads", ncores,
        extra.options
      ),
      verbose = verbose
    )

    rel_df <- bigreadr::fread2(paste0(prefix.out, ".kin0"), header = TRUE)
    names(rel_df) <- sub("^#(.*)$", "\\1", names(rel_df))
    rel_df

  }
}

################################################################################

#' Imputation
#'
#' Imputation using **Beagle** version 4.
#'
#' Downloads and more information can be found at the following websites
#' - [PLINK](http://www.cog-genomics.org/plink2),
#' - [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html).
#'
#' @param beagle.path Path to the executable of Beagle v4+.
#' @inheritParams snp_plinkRmSamples
#' @param bedfile.out Path to the output bedfile. Default is created by
#'   appending `"_impute"` to `prefix.in` (`bedfile.in` without extension).
#' @param memory.max Max memory (in GB) to be used. It is internally rounded
#'   to be an integer. Default is `3`.
#' @param extra.options Other options to be passed to Beagle as a string. More
#'   options can be found at Beagle's website.
#' @param plink.options Other options to be passed to PLINK as a string. More
#'   options can be found at \url{http://www.cog-genomics.org/plink2/filter}.
#' @param ncores Number of cores used. Default doesn't use parallelism.
#'   You may use [nb_cores].
#'
#' @references Browning, Brian L., and Sharon R. Browning.
#' "Genotype imputation with millions of reference samples."
#' The American Journal of Human Genetics 98.1 (2016): 116-126.
#'
#' @seealso [download_plink] [download_beagle]
#'
#' @return The path of the new bedfile.
#' @export
#'
snp_beagleImpute <- function(beagle.path,
                             plink.path,
                             bedfile.in,
                             bedfile.out = NULL,
                             memory.max = 3,
                             ncores = 1,
                             extra.options = "",
                             plink.options = "",
                             verbose = TRUE) {

  # get input and output right
  prefix.in <- sub_bed(bedfile.in)
  if (is.null(bedfile.out)) bedfile.out <- paste0(prefix.in, "_impute.bed")
  assert_noexist(bedfile.out)
  prefix.out <- sub_bed(bedfile.out)

  # get temporary files
  tmpfile1 <- tempfile()
  tmpfile2 <- tempfile()

  # Convert bed/bim/fam to vcf
  system_verbose(
    paste(
      plink.path,
      "--bfile", prefix.in,
      "--recode vcf bgz", # .vcf.gz extension
      "--out", tmpfile1,
      "--threads", ncores,
      plink.options
    ),
    verbose = verbose
  )
  vcf1 <- paste0(tmpfile1, ".vcf.gz")
  on.exit(file.remove(vcf1), add = TRUE)

  # Impute vcf with Beagle version 4
  system_verbose(
    paste(
      "java", sprintf("-Xmx%dg", round(memory.max)),
      "-jar", beagle.path,
      paste0("gt=", vcf1),
      paste0("out=", tmpfile2),
      paste0("nthreads=", ncores),
      extra.options
    ),
    verbose = verbose
  )
  vcf2 <- paste0(tmpfile2, ".vcf.gz")
  on.exit(file.remove(vcf2), add = TRUE)

  # Convert back vcf to bed/bim/fam
  system_verbose(
    paste(
      plink.path,
      "--vcf", vcf2,
      "--out", prefix.out,
      plink.options
    ),
    verbose = verbose
  )

  # return path to new bedfile
  bedfile.out
}

################################################################################
