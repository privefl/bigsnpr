% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autoSVD.R
\name{snp_autoSVD}
\alias{snp_autoSVD}
\alias{bed_autoSVD}
\title{Truncated SVD while limiting LD}
\usage{
snp_autoSVD(
  G,
  infos.chr,
  infos.pos = NULL,
  ind.row = rows_along(G),
  ind.col = cols_along(G),
  fun.scaling = snp_scaleBinom(),
  thr.r2 = 0.2,
  size = 100/thr.r2,
  k = 10,
  roll.size = 50,
  int.min.size = 20,
  alpha.tukey = 0.05,
  min.mac = 10,
  min.maf = 0.02,
  max.iter = 5,
  is.size.in.bp = NULL,
  ncores = 1,
  verbose = TRUE
)

bed_autoSVD(
  obj.bed,
  ind.row = rows_along(obj.bed),
  ind.col = cols_along(obj.bed),
  fun.scaling = bed_scaleBinom,
  thr.r2 = 0.2,
  size = 100/thr.r2,
  k = 10,
  roll.size = 50,
  int.min.size = 20,
  alpha.tukey = 0.05,
  min.mac = 10,
  min.maf = 0.02,
  max.iter = 5,
  ncores = 1,
  verbose = TRUE
)
}
\arguments{
\item{G}{A \link[bigstatsr:FBM.code256-class]{FBM.code256}
(typically \verb{<bigSNP>$genotypes}).\cr
\strong{You shouldn't have missing values.} Also, remember to do quality control,
e.g. some algorithms in this package won't work if you use SNPs with 0 MAF.}

\item{infos.chr}{Vector of integers specifying each SNP's chromosome.\cr
Typically \verb{<bigSNP>$map$chromosome}.}

\item{infos.pos}{Vector of integers specifying the physical position
on a chromosome (in base pairs) of each SNP.\cr
Typically \verb{<bigSNP>$map$physical.pos}.}

\item{ind.row}{An optional vector of the row indices (individuals) that
are used. If not specified, all rows are used.\cr
\strong{Don't use negative indices.}}

\item{ind.col}{An optional vector of the column indices (SNPs) that are used.
If not specified, all columns are used.\cr
\strong{Don't use negative indices.}}

\item{fun.scaling}{A function with parameters \code{X} (or \code{obj.bed}), \code{ind.row} and
\code{ind.col}, and that returns a data.frame with \verb{$center} and \verb{$scale} for the
columns corresponding to \code{ind.col}, to scale each of their elements such as followed:
\deqn{\frac{X_{i,j} - center_j}{scale_j}.} Default uses binomial scaling.
You can also provide your own \code{center} and \code{scale} by using \code{\link[bigstatsr:as_scaling_fun]{bigstatsr::as_scaling_fun()}}.}

\item{thr.r2}{Threshold over the squared correlation between two variants.
Default is \code{0.2}. Use \code{NA} if you want to skip the clumping step.}

\item{size}{For one SNP, window size around this SNP to compute correlations.
Default is \code{100 / thr.r2} for clumping (0.2 -> 500; 0.1 -> 1000; 0.5 -> 200).
If not providing \code{infos.pos} (\code{NULL}, the default), this is a window in
number of SNPs, otherwise it is a window in kb (physical distance).
I recommend that you provide the positions if available.}

\item{k}{Number of singular vectors/values to compute. Default is \code{10}.
\strong{This algorithm should be used to compute a few singular vectors/values.}}

\item{roll.size}{Radius of rolling windows to smooth log-p-values.
Default is \code{50}.}

\item{int.min.size}{Minimum number of consecutive outlier variants
in order to be reported as long-range LD region. Default is \code{20}.}

\item{alpha.tukey}{Default is \code{0.1}. The type-I error rate in outlier
detection (that is further corrected for multiple testing).}

\item{min.mac}{Minimum minor allele count (MAC) for variants to be included.
Default is \code{10}. Can actually be higher because of \code{min.maf}.}

\item{min.maf}{Minimum minor allele frequency (MAF) for variants to be included.
Default is \code{0.02}. Can actually be higher because of \code{min.mac}.}

\item{max.iter}{Maximum number of iterations of outlier detection.
Default is \code{5}.}

\item{is.size.in.bp}{Deprecated.}

\item{ncores}{Number of cores used. Default doesn't use parallelism.
You may use \code{\link[bigstatsr:reexports]{bigstatsr::nb_cores()}}.}

\item{verbose}{Output some information on the iterations? Default is \code{TRUE}.}

\item{obj.bed}{Object of type \link{bed}, which is the mapping of some bed file.
Use \code{obj.bed <- bed(bedfile)} to get this object.}
}
\value{
A named list (an S3 class "big_SVD") of
\itemize{
\item \code{d}, the singular values,
\item \code{u}, the left singular vectors,
\item \code{v}, the right singular vectors,
\item \code{niter}, the number of the iteration of the algorithm,
\item \code{nops}, number of Matrix-Vector multiplications used,
\item \code{center}, the centering vector,
\item \code{scale}, the scaling vector.
}

Note that to obtain the Principal Components, you must use
\link[bigstatsr:predict.big_SVD]{predict} on the result. See examples.
}
\description{
Fast truncated SVD with initial pruning and that iteratively removes
long-range LD regions. Some variants are removing due to the initial clumping,
then more and more variants are removed at each iteration. You can access the
indices of the remaining variants with \verb{attr(*, "subset")}. If some of the
variants removed are contiguous, the regions are reported in \verb{attr(*, "lrldr")}.
}
\details{
If you don't have any information about variants, you can try using
\itemize{
\item \code{infos.chr = rep(1, ncol(G))},
\item \code{size = ncol(G)} (if variants are not sorted),
\item \code{roll.size = 0} (if variants are not sorted).
}
}
\examples{
ex <- snp_attachExtdata()
G <- ex$genotypes

obj.svd <- snp_autoSVD(G,
                       infos.chr = ex$map$chromosome,
                       infos.pos = ex$map$physical.position)

str(obj.svd)

}
