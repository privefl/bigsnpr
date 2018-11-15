
#include <fstream>      // std::ifstream
#include <bigstatsr/utils.h>
#include <zlib.h>  // Need "PKG_LIBS = -lz"?
using namespace Rcpp;

int read_int(std::istream * ptr_is, std::streamsize n_byte = 4) {
  // https://stackoverflow.com/a/32066210/6103040
  int N;
  ptr_is->read(reinterpret_cast<char *>(&N), n_byte);
  return N;
}

std::string read_string(std::ifstream * ptr_is, std::streamsize n_byte = 2) {
  int len = read_int(ptr_is, n_byte);
  // https://stackoverflow.com/a/38623543/6103040
  char buffer[len + 1];
  ptr_is->read(buffer, len);
  buffer[len] = '\0';
  return std::string(buffer, len);
}

struct membuf : std::streambuf
{
  membuf(char* begin, char* end) {
    this->setg(begin, begin, end);
  }
};

// [[Rcpp::export]]
void read_bgen(std::string filename, int N = 487409) {

  std::ifstream is(filename.c_str(), std::ifstream::binary);
  if (!is) Rcpp::stop("Error while opening '%s'.", filename);

  // std::istream::sentry se(is, true);
  // std::streambuf* sb = is.rdbuf();

  int offset = read_int(&is);
  Rcout << offset << std::endl;

  // Get past the header block
  is.seekg(offset, std::ios_base::cur);

  std::string id = read_string(&is);
  Rcout << id << std::endl;

  std::string rsid = read_string(&is);
  Rcout << rsid << std::endl;

  std::string chr = read_string(&is);
  Rcout << chr << std::endl;

  int pos = read_int(&is);
  Rcout << pos << std::endl;

  int K = read_int(&is, 2);
  myassert(K == 2, "Only 2 alleles allowed.");

  std::string a1 = read_string(&is, 4);
  Rcout << a1 << std::endl;
  std::string a2 = read_string(&is, 4);
  Rcout << a2 << std::endl;

  int C = read_int(&is);
  Rcout << C << std::endl;
  int D = read_int(&is);
  myassert(D == (10 + 3 * N), "Wrong size of decompressed data.");

  unsigned char buffer_in[C];
  unsigned char buffer_out[D];
  is.read((char *)buffer_in, C);

  // https://gist.github.com/arq5x/5315739
  // zlib struct
  z_stream infstream;
  infstream.zalloc = Z_NULL;
  infstream.zfree = Z_NULL;
  infstream.opaque = Z_NULL;
  infstream.avail_in = C;            // size of input
  infstream.next_in = buffer_in;     // input char array
  infstream.avail_out = D;           // size of output
  infstream.next_out = buffer_out;   // output char array
  // the actual DE-compression work.
  myassert(inflateInit(&infstream) == Z_OK, "Problem when decompressing.");
  myassert(inflate(&infstream, Z_NO_FLUSH) != Z_STREAM_ERROR,
           "Problem when decompressing.");
  inflateEnd(&infstream);

  membuf sbuf((char *)buffer_out, (char *)buffer_out + sizeof(buffer_out));
  std::istream in(&sbuf);
  myassert(read_int(&in) == N, "Problem with N.");
  myassert(read_int(&in, 2) == K, "Problem with K.");
  myassert(read_int(&in, 1) == K, "Problem with P_min.");
  myassert(read_int(&in, 1) == K, "Problem with P_max.");
  in.seekg(N, std::ios_base::cur); // Skip ploidy
  myassert(read_int(&in, 1) == 0, "Problem with 'phased'.");

  is.close();
}

/*** R
file <- "tmp-data/first_bytes.bgen"
if (!file.exists(file)) file <- paste0("../", file)
read_bgen(file)
*/
