################################################################################

#' Replace extension '.bed'
#'
#' @param path String with extension '.bed'.
#' @param replacement Replacement of '.bed'. Default replaces by nothing.
#'   Can be useful to replace e.g. by '.bim' or '.fam'.
#' @param stop_if_not_ext If `replacement != ""`, whether to error if
#'   replacement is not an extension (starting with a '.').
#'
#' @return String with extension '.bed' replaced by `replacement`.
#' @export
#'
#' @examples
#' path <- "toto.bed"
#' sub_bed(path)
#' sub_bed(path, ".bim")
#' sub_bed(path, ".fam")
#' sub_bed(path, "_QC", stop_if_not_ext = FALSE)
sub_bed <- function(path, replacement = "", stop_if_not_ext = TRUE) {
  pattern <- "\\.bed$"
  if (!grepl(pattern, path))
    stop2("Path '%s' must have 'bed' extension.", path)
  if (stop_if_not_ext && nchar(replacement) > 0 && substr(replacement, 1, 1) != ".")
    stop2("Replacement must be an extension starting with '.' if provided.")
  sub(pattern, replacement, path)
}

################################################################################

#' Class bed
#'
#' A reference class for storing a pointer to a mapped version of a bed file.
#'
#' @details
#' A `bed` object has many field:
#'   - `$address`: address of the external pointer containing the underlying
#'     C++ object, to be used as a `XPtr<bed>` in C++ code
#'   - `$extptr`: use `$address` instead
#'   - `$bedfile`: path to the bed file
#'   - `$bimfile`: path to the corresponding bim file
#'   - `$famfile`: path to the corresponding fam file
#'   - `$nrow`: number of samples in the bed file
#'   - `$ncol`: number of variants in the bed file
#'   - `$map`: data frame read from `$bimfile`
#'   - `$fam`: data frame read from `$famfile`
#'   - `$.map`: use `$map` instead
#'   - `$.fam`: use `$fam` instead
#'
#' @examples
#' bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
#' obj.bed <- bed(bedfile)
#'
#' @exportClass bed
#' @importFrom methods new
#'
bed_RC <- methods::setRefClass(

  "bed",

  fields = list(
    bedfile = "character",
    extptr  = "externalptr",
    .fam    = "data.frame",
    .map    = "data.frame",

    #### Active bindings
    bimfile = function() sub_bed(.self$bedfile, ".bim"),
    famfile = function() sub_bed(.self$bedfile, ".fam"),

    fam = function() {
      if (base::ncol(.self$.fam) == 0) {
        .self$.fam <- bigreadr::fread2(.self$famfile, col.names = NAMES.FAM)
      }
      .self$.fam
    },
    map = function() {
      if (base::ncol(.self$.map) == 0) {
        .self$.map <- bigreadr::fread2(.self$bimfile, col.names = NAMES.MAP)
      }
      .self$.map
    },

    nrow = function() base::nrow(.self$fam),
    ncol = function() base::nrow(.self$map),

    # infos.chr = function() bigreadr::fread2(.self$bimfile, select = 1)[[1]],
    # infos.pos = function() bigreadr::fread2(.self$bimfile, select = 4)[[1]],

    address = function() {
      if (identical(.self$extptr, new("externalptr"))) { # nil
        .self$extptr <- bedXPtr(.self$bedfile, .self$nrow, .self$ncol)
      }
      .self$extptr
    }
  ),

  methods = list(
    initialize = function(bedfile) {

      .self$bedfile <- path.expand(bedfile)

      # Check if all three files exist
      sapply(c(.self$bedfile, .self$bimfile, .self$famfile), assert_exist)

      # Connect once
      .self$address

      .self
    },

    show = function() {
      cat(sprintf(
        "A 'bed' object with %s samples and %s variants.\n",
        .self$nrow, .self$ncol))
      invisible(.self)
    }
  )
)

################################################################################

#' Wrapper constructor for class `bed`.
#'
#' @inheritParams snp_readBed
#'
#' @rdname bed-class
#'
#' @export
#'
bed <- function(bedfile) new(Class = "bed", bedfile = bedfile)

################################################################################

#' Methods for the bed class
#'
#' @param x Object of type `bed`.
#'
#' @name bed-methods
NULL

################################################################################

#' Dimension methods for class `bed`.
#' Methods `nrow()` and `ncol()` are automatically defined with `dim()`.
#'
#' @rdname bed-methods
#' @export
setMethod("dim",    signature(x = "bed"), function(x) c(x$nrow, x$ncol))

#' @rdname bed-methods
#' @export
setMethod("length", signature(x = "bed"), function(x) prod(dim(x) + 0))

################################################################################
