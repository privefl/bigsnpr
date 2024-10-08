% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binom-scaling.R
\name{bed_counts}
\alias{bed_counts}
\title{Counts}
\usage{
bed_counts(
  obj.bed,
  ind.row = rows_along(obj.bed),
  ind.col = cols_along(obj.bed),
  byrow = FALSE,
  ncores = 1
)
}
\arguments{
\item{obj.bed}{Object of type \link{bed}, which is the mapping of some bed file.
Use \code{obj.bed <- bed(bedfile)} to get this object.}

\item{ind.row}{An optional vector of the row indices (individuals) that
are used. If not specified, all rows are used.\cr
\strong{Don't use negative indices.}}

\item{ind.col}{An optional vector of the column indices (SNPs) that are used.
If not specified, all columns are used.\cr
\strong{Don't use negative indices.}}

\item{byrow}{Whether to count by individual rather than by variant?
Default is \code{FALSE} (count by variant).}

\item{ncores}{Number of cores used. Default doesn't use parallelism.
You may use \code{\link[bigstatsr:reexports]{bigstatsr::nb_cores()}}.}
}
\value{
A matrix of with 4 rows and \code{length(ind.col)} columns.
}
\description{
Counts the number of 0s, 1s, 2s and NAs by variants in the bed file.
}
\examples{
bedfile <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
obj.bed <- bed(bedfile)

bed_counts(obj.bed, ind.col = 1:5)

bed_counts(obj.bed, ind.row = 1:5, byrow = TRUE)

}
