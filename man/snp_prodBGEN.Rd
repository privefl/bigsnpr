% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prod-bgen.R
\name{snp_prodBGEN}
\alias{snp_prodBGEN}
\title{BGEN matrix product}
\usage{
snp_prodBGEN(
  bgenfiles,
  beta,
  list_snp_id,
  ind_row = NULL,
  bgi_dir = dirname(bgenfiles),
  read_as = c("dosage", "random"),
  block_size = 1000,
  ncores = 1
)
}
\arguments{
\item{bgenfiles}{Character vector of paths to files with extension ".bgen".
The corresponding ".bgen.bgi" index files must exist.}

\item{beta}{A matrix (or a vector), with rows corresponding to \code{list_snp_id}.}

\item{list_snp_id}{List of character vectors of SNP IDs to read, with one
vector per BGEN file. Each SNP ID should be in the form
\code{"<chr>_<pos>_<a1>_<a2>"} (e.g. \code{"1_88169_C_T"} or \code{"01_88169_C_T"}).
If you have one BGEN file only, just wrap your vector of IDs with \code{list()}.
\strong{This function assumes that these IDs are uniquely identifying variants.}}

\item{ind_row}{An optional vector of the row indices (individuals) that
are used. If not specified, all rows are used. \strong{Don't use negative indices.}
You can access the sample IDs corresponding to the genotypes from the \emph{.sample}
file, and use e.g. \code{match()} to get indices corresponding to the ones you want.}

\item{bgi_dir}{Directory of index files. Default is the same as \code{bgenfiles}.}

\item{read_as}{How to read BGEN probabilities? Currently implemented:
\itemize{
\item as dosages (rounded to two decimal places), the default,
\item as hard calls, randomly sampled based on those probabilities
(similar to PLINK option '\verb{--hard-call-threshold random}').
}}

\item{block_size}{Maximum size of temporary blocks (in number of variants).
Default is \code{1000}.}

\item{ncores}{Number of cores used. Default doesn't use parallelism.
You may use \code{\link[bigstatsr:reexports]{bigstatsr::nb_cores()}}.}
}
\value{
The product \code{bgen_data[ind_row, 'list_snp_id'] \%*\% beta}.
}
\description{
Compute a matrix product between BGEN files and a matrix. This removes the
need to read an intermediate FBM object with \code{\link[=snp_readBGEN]{snp_readBGEN()}} to compute the
product. Moreover, when using dosages, they are not rounded to two decimal
places anymore.
}
\seealso{
\code{\link[=snp_readBGEN]{snp_readBGEN()}}
}
