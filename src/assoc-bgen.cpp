/******************************************************************************/

#include <bigstatsr/BMAcc.h>
#include <zlib.h>
#include "fstream.h"

using namespace Rcpp;

/******************************************************************************/

void r2_variant(std::ifstream * ptr_stream,
                const LogicalVector& use_ind,
                const NumericVector& decode,
                const NumericVector& y,
                double ySum,
                double * ptr_r2) {

  std::string id   = read_string(ptr_stream);
  std::string rsid = read_string(ptr_stream);
  std::string chr  = read_string(ptr_stream);
  int pos = read_int(ptr_stream);
  myassert(pos > 0, "Positions should be positive.");
  int K   = read_int(ptr_stream, 2);
  myassert(K == 2, "Only 2 alleles allowed.");
  std::string a1 = read_string(ptr_stream, 4);
  std::string a2 = read_string(ptr_stream, 4);

  int C = read_int(ptr_stream) - 4;
  int D = read_int(ptr_stream);

  // decompress variant data
  unsigned char *buffer_in = new unsigned char[C];
  ptr_stream->read((char *)buffer_in, C);
  unsigned char *buffer_out = new unsigned char[D];
  // zlib struct (https://gist.github.com/arq5x/5315739)
  z_stream infstream;
  infstream.zalloc = Z_NULL;
  infstream.zfree  = Z_NULL;
  infstream.opaque = Z_NULL;
  infstream.avail_in  = C;             // size of input
  infstream.next_in   = buffer_in;     // input char array
  infstream.avail_out = D;             // size of output
  infstream.next_out  = buffer_out;    // output char array
  // the actual DE-compression work.
  myassert(inflateInit(&infstream) == Z_OK, "Problem when decompressing.");
  myassert(inflate(&infstream, Z_NO_FLUSH) != Z_STREAM_ERROR,
           "Problem when decompressing.");
  inflateEnd(&infstream);

  // read decompress "probabilities" as dosages
  // and compute squared correlation with y
  int i = 0, i2 = 0, i_pld = 8, N = (D - 10) / 3, i_prblt = 10 + N, i_code;
  double x, num, deno_x, xSum = 0, xxSum = 0, xySum = 0;
  for (; i_prblt < D; i++, i_pld++, i_prblt += 2) {
    if (use_ind[i]) {
      if (buffer_out[i_pld] >= 0x80) {
        return; // at least one missing value
      } else {
        i_code = 2 * buffer_out[i_prblt] + buffer_out[i_prblt + 1];
        x = decode[i_code];
        xSum  += x;
        xxSum += x * x;
        xySum += x * y[i2++];
      }
    }
  }
  int n = y.size();
  num = xySum - xSum * ySum / n;
  deno_x = xxSum - xSum * xSum / n;
  *ptr_r2 = num * num / deno_x;  // divide by deno_y later

  delete[] buffer_in;
  delete[] buffer_out;
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector r2_bgen(std::string filename,
                      NumericVector offsets,
                      LogicalVector use_ind,
                      NumericVector decode,
                      NumericVector y) {

  int K = offsets.size();
  NumericVector r2(K, NA_REAL);
  double y_sum = Rcpp::sum(y);

  // connection to BGEN file
  std::ifstream stream(filename.c_str(), std::ifstream::binary);
  if (!stream) Rcpp::stop("Error while opening '%s'.", filename);

  // read variants one by one
  for (int k = 0; k < K; k++) {
    stream.seekg(offsets[k]);
    r2_variant(&stream, use_ind, decode, y, y_sum, &(r2[k]));
  }

  stream.close();

  return r2;
}

/******************************************************************************/
