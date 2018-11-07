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
                  boost::unordered_map<std::string, int>& mymap,
                  BMAcc<unsigned char>& macc,
                  CharacterVector& ID,
                  CharacterVector& RSID,
                  CharacterVector& CHR,
                  NumericVector  & POS,
                  CharacterVector& A1,
                  CharacterVector& A2) {

  std::string id = read_string(ptr_stream);
  boost::unordered_map<std::string, int>::iterator got = mymap.find(id);

  if (got == mymap.end()) {

    // we don't want this variant -> go to next variant
    ptr_stream->seekg(18, std::ios_base::cur);
    int C = read_int(ptr_stream);
    ptr_stream->seekg(C, std::ios_base::cur);

  } else {

    Rcout << "Found one!!" << std::endl;

    // we want this variant -> get info
    std::string rsid = read_string(ptr_stream);
    std::string chr  = read_string(ptr_stream);
    int pos = read_int(ptr_stream);
    int K   = read_int(ptr_stream, 2);
    myassert(K == 2, "Only 2 alleles allowed.");
    std::string a1 = read_string(ptr_stream, 4);
    std::string a2 = read_string(ptr_stream, 4);

    // store variant info
    int j = got->second;
    ID[j]   = id;
    RSID[j] = rsid;
    CHR[j]  = chr;
    POS[j]  = pos;
    A1[j]   = a1;
    A2[j]   = a2;

    // decompress variant data
    int C = read_int(ptr_stream) - 4;
    int D = read_int(ptr_stream);
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
    double x, coeff = 100 / 255;
    int N = (D - 10) / 3;
    for (int i = 0, i2 = 10 + N; i2 < N; i++, i2 += 2) { // Skip infos + ploidy
      x = 2 * buffer_out[i2] + buffer_out[i2 + 1];
      macc(i, j) = 207 - std::round(coeff * x);
    }

    // no longer need to search for this variant
    mymap.erase(got);

  }
}

/******************************************************************************/

void read_file(std::string filename,
               boost::unordered_map<std::string, int>& mymap,
               BMAcc<unsigned char>& macc,
               CharacterVector& ID,
               CharacterVector& RSID,
               CharacterVector& CHR,
               NumericVector  & POS,
               CharacterVector& A1,
               CharacterVector& A2) {

  std::ifstream stream(filename.c_str(), std::ifstream::binary);
  if (!stream) Rcpp::stop("Error while opening '%s'.", filename);

  // get past the header block -> go to first variant
  int offset = read_int(&stream);
  read_int(&stream);
  int n_var = read_int(&stream);  // number of variants in the file
  stream.seekg(offset - 8, std::ios_base::cur);

  for (int j = 0; j < n_var; j++) {
    if (j % 100000 == 0) {
      Rcout << j << std::endl;
    }
    read_variant(&stream, mymap, macc, ID, RSID, CHR, POS, A1, A2);
  }

  stream.close();
}

/******************************************************************************/

// [[Rcpp::export]]
DataFrame readbgen(const CharacterVector& filenames,
                   const CharacterVector& snp_id,
                   Environment BM) {

  // FBM accessor
  XPtr<FBM> xpBM = BM["address"];
  BMAcc<unsigned char> macc(xpBM);

  // output variant info
  int M = snp_id.size();
  CharacterVector ID(M, NA_STRING), RSID(M, NA_STRING), CHR(M, NA_STRING);
  NumericVector POS(M, NA_REAL);
  CharacterVector A1(M, NA_STRING), A2(M, NA_STRING);

  // use boost::unordered_map to speed up search
  boost::unordered_map<std::string, int> mymap;
  mymap.reserve(M);
  for (int j = 0; j < M; j++) {
    mymap.insert(std::make_pair(as<std::string>(snp_id[j]), j));
  }

  // read BGEN files one by one
  for (int k = 0; k < filenames.size(); k++) {
    read_file(as<std::string>(filenames[k]), mymap, macc,
              ID, RSID, CHR, POS, A1, A2);
  }

  // warn if some variants have not been found
  if (mymap.size() > 0) {
    Rcpp::warning("%d variants have not been matched.", mymap.size());
    int N = macc.nrow();
    boost::unordered_map<std::string, int>::iterator it = mymap.begin(),
      end = mymap.end();
    while (it != end) {
      int j = it->second;
      for (int i = 0; i < N; i++) macc(i, j) = 3;  // set as missing
      it++;
    }
  }

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
