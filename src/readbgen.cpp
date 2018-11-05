/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <boost/unordered_map.hpp>
#include <fstream>
#include <zlib.h>

using namespace Rcpp;

/******************************************************************************/

int read_int(std::istream * ptr_stream, std::streamsize n_byte = 4) {
  // https://stackoverflow.com/a/32066210/6103040
  int N;
  ptr_stream->read(reinterpret_cast<char *>(&N), n_byte);
  return N;
}

std::string read_string(std::ifstream * ptr_stream, std::streamsize n_byte = 2) {
  int len = read_int(ptr_stream, n_byte);
  // https://stackoverflow.com/a/38623543/6103040
  char buffer[len + 1];
  ptr_stream->read(buffer, len);
  buffer[len] = '\0';
  return std::string(buffer, len);
}

/******************************************************************************/

void read_variant(std::ifstream * ptr_stream,
                  boost::unordered_map<std::string, int>& mymap,
                  CharacterVector& ID,
                  CharacterVector& RSID,
                  CharacterVector& CHR,
                  NumericVector&   POS,
                  CharacterVector& A1,
                  CharacterVector& A2) {

  boost::unordered_map<std::string, int>::iterator got;

  std::string id = read_string(ptr_stream);
  std::string rsid = read_string(ptr_stream);
  std::string chr = read_string(ptr_stream);
  int pos = read_int(ptr_stream);
  int K = read_int(ptr_stream, 2);
  myassert(K == 2, "Only 2 alleles allowed.");
  std::string a1 = read_string(ptr_stream, 4);
  std::string a2 = read_string(ptr_stream, 4);

  int C = read_int(ptr_stream) - 4;
  int D = read_int(ptr_stream);

  got = mymap.find(id);
  if (got == mymap.end()) {

    // we don't want this variant -> go to next variant
    ptr_stream->seekg(C, std::ios_base::cur);

  } else {

    // we want this variant
    int j = got->second;

    // store variant info
    ID[j]   = id;
    RSID[j] = rsid;
    CHR[j]  = chr;
    POS[j]  = pos;
    A1[j]   = a1;
    A2[j]   = a2;

    // decompress variant data
    unsigned char buffer_in[C];
    ptr_stream->read((char *)buffer_in, C);
    unsigned char buffer_out[D];
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

    // read decompress data
    double x, coeff = 100 / 255;
    int N = (D - 10) / 3;
    for (int i = 0, i2 = 10 + N; i2 < N; i++, i2 += 2) { // Skip infos + ploidy
      x = 2 * buffer_out[i2] + buffer_out[i2 + 1];
      macc(i, j) = 207 - std::round(coeff * x);
    }

  }
}

/******************************************************************************/

void read_file(std::string filename) {

  std::ifstream stream(filename.c_str(), std::ifstream::binary);
  if (!stream) Rcpp::stop("Error while opening '%s'.", filename);

  // Get past the header block -> go to first variant
  int offset = read_int(&stream);
  read_int(&stream);
  int n_var = read_int(&stream);  // number of variants in the file
  stream.seekg(offset - 8, std::ios_base::cur);

  for (int j = 0; j < n_var; j++) {
    read_variant(&stream);
  }

  stream.close();
}

/******************************************************************************/

// [[Rcpp::export]]
DataFrame readbgen(const CharacterVector& filenames,
                   const CharacterVector& snp_id,
                   Environment BM) {

  XPtr<FBM> xpBM = BM["address"];
  unsigned char* ptr = static_cast<unsigned char*>(xpBM->matrix());

  int M = snp_id.size();
  CharacterVector ID(M, NA_STRING), RSID(M, NA_STRING), CHR(M, NA_STRING);
  NumericVector POS(M, NA_REAL);
  CharacterVector A1(M, NA_STRING), A2(M, NA_STRING);

  for (int k = 0; k < filenames.size(); k++) {
    read_file(filenames[k], ID, RSID, CHR, POS, A1, A2);
  }

  // Warn if some SNPs have not been found

  return DataFrame::create(
    Named("chromosome") = CHR,
    Named("marker.ID") = ID,
    Named("rsID") = RSID,
    Named("physical.pos") = POS,
    Named("allele1") = A1,
    Named("allele2") = A2
  );
}

/******************************************************************************/
