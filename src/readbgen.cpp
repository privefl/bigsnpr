/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <boost/unordered_map.hpp>
#include <fstream>
#include <zlib.h>

using namespace Rcpp;

/******************************************************************************/

// https://stackoverflow.com/a/32066210/6103040
inline int read_uint32(std::istream * ptr_stream) {
  uint32_t N;
  ptr_stream->read(reinterpret_cast<char *>(&N), 4);
  return N;
}

inline int read_uint16(std::istream * ptr_stream) {
  uint16_t N;
  ptr_stream->read(reinterpret_cast<char *>(&N), 2);
  return N;
}

/******************************************************************************/

inline int read_int(std::istream * ptr_stream, std::streamsize n_byte = 4) {
  if (n_byte == 4) {
    return read_uint32(ptr_stream);
  } else if (n_byte == 2) {
    return read_uint16(ptr_stream);
  } else {
    Rcpp::stop("Not supported.");
  }
}

inline std::string read_string(std::ifstream * ptr_stream,
                               std::streamsize n_byte = 2) {
  int len = read_int(ptr_stream, n_byte);
  // https://stackoverflow.com/a/38623543/6103040
  char buffer[len + 1];
  ptr_stream->read(buffer, len);
  buffer[len] = '\0';
  return std::string(buffer, len);
}

/******************************************************************************/

void read_variant(std::ifstream * ptr_stream,
                  unsigned char * ptr_mat,
                  const RawVector& decode) {

  std::string id   = read_string(ptr_stream);
  std::string rsid = read_string(ptr_stream);
  Rcout << rsid << std::endl;
  std::string chr  = read_string(ptr_stream);
  int pos = read_int(ptr_stream);
  int K   = read_int(ptr_stream, 2);
  // myassert(K == 2, "Only 2 alleles allowed.");
  std::string a1 = read_string(ptr_stream, 4);
  std::string a2 = read_string(ptr_stream, 4);

  int C = read_int(ptr_stream) - 4;
  int D = read_int(ptr_stream);

  // decompress variant data
  unsigned char buffer_in[C];
  ptr_stream->read((char *)buffer_in, C);
  unsigned char buffer_out[D];
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

  // read decompress probabilities and store them as rounded dosages
  int x;
  int N = (D - 10) / 3;
  for (int i = 0, i2 = 10 + N; i2 < D; i++, i2 += 2) { // Skip infos + ploidy
    x = 2 * buffer_out[i2] + buffer_out[i2 + 1];
    ptr_mat[i] = decode[x];
  }
}

/******************************************************************************/

// [[Rcpp::export]]
void read_bgen(std::string filename,
               NumericVector offsets,
               Environment BM,
               IntegerVector ind_col,
               RawVector decode) {

  XPtr<FBM> xpBM = BM["address"];
  unsigned char* ptr_mat = static_cast<unsigned char*>(xpBM->matrix());

  std::size_t j, n = xpBM->nrow();
  int K = offsets.size();
  myassert(ind_col.size() == K, ERROR_DIM);

  // connection to BGEN file
  std::ifstream stream(filename.c_str(), std::ifstream::binary);
  if (!stream) Rcpp::stop("Error while opening '%s'.", filename);

  // read variants one by one
  for (int k = 0; k < K; k++) {
    stream.seekg(offsets[k]);
    j = ind_col[k] - 1;
    read_variant(&stream, ptr_mat + n * j, decode);
  }

  stream.close();
}

/******************************************************************************/
