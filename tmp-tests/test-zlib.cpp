#include <Rcpp.h>
using namespace Rcpp;

#include <zlib.h>  // Need "PKG_LIBS = -lz"?
#include <Rdefines.h>
#include <Rinternals.h>
#include "R_ext/Memory.h"
#include "R_ext/Utils.h"

static voidpf R_zlib_alloc(voidpf ptr, uInt items, uInt size) {
  return(R_alloc(items, size));
}

static void R_zlib_free(voidpf ptr, voidpf addr) {}

//' Inflate a raw vector that was compressed
//'
//' @param r_source raw vector of compressed data.
//' @param format A character scalar indicating the type of compression that was applied.
//' This must be one of "gzip", "zlib" or "raw".
//' @param r_guess_size your best guess as to the size of uncompressed data. Not ideal,
//'        and this won't be necessary in future releases. Reember, this is a direct
//'        port from the defunct \code{Rcompression} package.\cr
//'        \cr
//'        Aim high as you'll only get back the actual number of bytes
//'        in the uncompressed data.
//' @return Raw vector of expanded data.
//' @export
//' @examples
//' raw <- charToRaw("The quick brown fox jumps over the lazy dog.")
//' compressed <- mem_compress(raw, format = "gzip")
//' decompressed <- mem_inflate(compressed, format = "gzip", 1000)
// [[Rcpp::export]]
SEXP mem_inflate(SEXP r_source, String format, SEXP r_guess_size) {

  z_stream stream;
  int err, guess_size = REAL(r_guess_size)[0], windowBits;
  unsigned char *ans = (unsigned char *)R_alloc(guess_size, sizeof(unsigned char));
  SEXP r_ans;

  if (format == "gzip"){
    windowBits = MAX_WBITS+16;
  } else if (format == "zlib"){
    windowBits = MAX_WBITS;
  } else if (format == "raw"){
    windowBits = -MAX_WBITS;
  } else {
    Rcpp::stop("Invalid format argument");
  }

  stream.next_in = RAW(r_source);
  stream.avail_in = GET_LENGTH(r_source);
  stream.next_out = (unsigned char *)ans;
  stream.avail_out = guess_size;

  stream.zalloc = R_zlib_alloc;
  stream.zfree = R_zlib_free;
  stream.opaque = NULL;

  err = inflateInit2(&stream, windowBits);
  if(err != Z_OK) {
    Rcpp::stop("cannot establish the uncompress/inflate stream on this data");
  }

  err = inflate(&stream, Z_FINISH);

  if (err < 0) {
    inflateEnd(&stream);
    Rcpp::stop("Failed to uncompress the raw data");
  }

  inflateEnd(&stream);

  r_ans = Rf_allocVector(RAWSXP, stream.total_out);
  memcpy(RAW(r_ans), ans, stream.total_out);

  return(r_ans);

}

// [[Rcpp::export]]
SEXP mem_inflate2(RawVector r_source, int guess_size) {

  z_stream stream;
  int err;
  unsigned char *ans = (unsigned char *)R_alloc(guess_size, sizeof(unsigned char));
  SEXP r_ans;

  stream.next_in = RAW(r_source);
  stream.avail_in = GET_LENGTH(r_source);
  stream.next_out = (unsigned char *)ans;
  stream.avail_out = guess_size;

  stream.zalloc = R_zlib_alloc;
  stream.zfree = R_zlib_free;
  stream.opaque = NULL;

  err = inflateInit2(&stream, MAX_WBITS);
  if(err != Z_OK) {
    Rcpp::stop("cannot establish the uncompress/inflate stream on this data");
  }

  err = inflate(&stream, Z_FINISH);

  if (err < 0) {
    inflateEnd(&stream);
    Rcpp::stop("Failed to uncompress the raw data");
  }

  inflateEnd(&stream);

  r_ans = Rf_allocVector(RAWSXP, stream.total_out);
  memcpy(RAW(r_ans), ans, stream.total_out);

  return(r_ans);
}

/*** R
N <- 1e3
data <- sample(as.raw(0:100), size = N, replace = TRUE)
compressed <- memCompress(data, type = "gzip")
mem_inflate(compressed, format = "zlib", r_guess_size = N)
mem_inflate2(compressed, guess_size = N)

microbenchmark::microbenchmark(
  mem_inflate(compressed, format = "zlib", r_guess_size = N),
  mem_inflate2(compressed, guess_size = N),
  times = 200
)
*/
