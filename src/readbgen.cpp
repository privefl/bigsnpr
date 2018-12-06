/******************************************************************************/

#include <bigstatsr/BMAcc.h>
#include <zlib.h>
#include "fstream.h"

using namespace Rcpp;

/******************************************************************************/

void read_variant(std::ifstream * ptr_stream,
                  unsigned char * ptr_mat,
                  const IntegerVector& ind_row,
                  const RawVector& decode,
                  std::vector<std::string>& ID) {

  std::string id   = read_string(ptr_stream);
  ID.push_back(id);
  std::string rsid = read_string(ptr_stream);
  // Rcout << rsid << std::endl;
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

  // read decompress "probabilities" and store them as rounded dosages
  int i, i_pld, i_prblt, i_G, x, N = (D - 10) / 3;
  for (i = 0, i_pld = 8, i_prblt = 10 + N; i_prblt < D; i++, i_pld++, i_prblt += 2) {
    i_G = ind_row[i];
    if (i_G >= 0) {
      if (buffer_out[i_pld] >= 0x80) {
        ptr_mat[i_G] = 208;  // missing
      } else {
        x = 2 * buffer_out[i_prblt] + buffer_out[i_prblt + 1];
        ptr_mat[i_G] = decode[x];
      }
    }
  }
}

/******************************************************************************/

// [[Rcpp::export]]
CharacterVector read_bgen(std::string filename,
                          NumericVector offsets,
                          Environment BM,
                          IntegerVector ind_row,
                          IntegerVector ind_col,
                          RawVector decode) {

  XPtr<FBM> xpBM = BM["address"];
  unsigned char* ptr_mat = static_cast<unsigned char*>(xpBM->matrix());

  std::size_t j, n = xpBM->nrow();
  int K = offsets.size();
  myassert(ind_col.size() == K, ERROR_DIM);
  std::vector<std::string> ID; ID.reserve(K);

  // connection to BGEN file
  std::ifstream stream(filename.c_str(), std::ifstream::binary);
  if (!stream) Rcpp::stop("Error while opening '%s'.", filename);

  // read variants one by one
  for (int k = 0; k < K; k++) {
    // if (k % 10000 == 0) Rcout << k << std::endl;
    stream.seekg(offsets[k]);
    j = ind_col[k] - 1;
    read_variant(&stream, ptr_mat + n * j, ind_row, decode, ID);
  }

  stream.close();

  return wrap(ID);
}

/******************************************************************************/
