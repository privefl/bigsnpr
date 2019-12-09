#include <Rcpp.h>
#include <fstream>
#include <zlib.h>

using namespace Rcpp;

inline void myassert(bool cond, const char *msg) {
  if (!cond) Rcpp::stop(msg);
}

// [[Rcpp::export]]
void read_compressed(std::string filename, size_t size = 64) {

  std::ifstream stream(filename.c_str(), std::ifstream::binary);
  if (!stream) Rcpp::stop("Error while opening '%s'.", filename);

  int C = size;
  int D = size * 10;

  // decompress variant data
  unsigned char *buffer_in = new unsigned char[C];
  stream.read((char *)buffer_in, C);
  unsigned char *buffer_out = new unsigned char[D];
  // zlib struct (https://gist.github.com/arq5x/5315739)
  z_stream infstream;
  infstream.zalloc = Z_NULL;
  infstream.zfree  = Z_NULL;
  infstream.opaque = Z_NULL;
  infstream.data_type = Z_TEXT;
  infstream.avail_in  = C;             // size of input
  infstream.next_in   = buffer_in;     // input char array
  infstream.avail_out = D;             // size of output
  infstream.next_out  = buffer_out;    // output char array
  // the actual DE-compression work.
  myassert(inflateInit(&infstream) == Z_OK, "Problem when decompressing.");
  myassert(inflate(&infstream, Z_NO_FLUSH) != Z_STREAM_ERROR,
           "Problem when decompressing.");
  inflateEnd(&infstream);

  Rcpp::Rcout << buffer_out << std::endl;

  delete[] buffer_in;
  delete[] buffer_out;
  stream.close();
}

/*** R
capture.output(read_compressed("tmp-data/test.txt.gz", 1000))
(.Last.value)
*/
