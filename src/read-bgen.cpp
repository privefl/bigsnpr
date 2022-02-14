/******************************************************************************/

#include <bigstatsr/BMAcc.h>
#include <zlib.h>
#include "fstream.h"

using namespace Rcpp;

/******************************************************************************/

inline unsigned char sample_from_prob_uchar(double p0, double p1) {
  double first = ::unif_rand() * 255 - p0;
  return((first < 0) ? 4 : ((first < p1) ? 5 : 6));
}

/******************************************************************************/

std::string read_variant(std::ifstream * ptr_stream,
                         unsigned char * ptr_mat,
                         const IntegerVector& ind_row,
                         const RawVector& decode,
                         bool dosage,
                         int N,
                         double * info,
                         double * freq) {

  std::string id   = read_string(ptr_stream);
  std::string rsid = read_string(ptr_stream);
  std::string chr  = read_string(ptr_stream);
  int pos = read_int(ptr_stream);
  int K   = read_int(ptr_stream, 2);
  myassert(pos > 0, "Positions should be positive.");
  myassert(K == 2,  "Only 2 alleles allowed.");
  std::string a1 = read_string(ptr_stream, 4);
  std::string a2 = read_string(ptr_stream, 4);

  int C = read_int(ptr_stream) - 4;
  int D = read_int(ptr_stream);
  myassert(D == (10 + 3 * N), "Probabilities should be stored using 8 bits.");

  // uncompress variant data
  unsigned char *buffer_in = new unsigned char[C];
  ptr_stream->read((char *)buffer_in, C);
  unsigned char *buffer_out = new unsigned char[D];
  uLongf D2 = D;
  myassert(uncompress(buffer_out, &D2, buffer_in, C) == Z_OK,
           "Problem when uncompressing.");

  // read uncompressed "probabilities" and store them as rounded dosages or hard calls
  int n = ind_row.size();
  int nona = n;
  double af = 0, num = 0;
  for (int i = 0; i < n; i++) {
    int i_G = ind_row[i];
    int i_pld = 8 + i_G;
    if (buffer_out[i_pld] >= 0x80) {
      ptr_mat[i] = 3;  // missing
      nona--;
    } else {
      int i_prblt = 10 + N + 2 * i_G;
      // probabilities * 255
      unsigned char p0 = buffer_out[i_prblt];
      unsigned char p1 = buffer_out[i_prblt + 1];
      double e_ij = 2 * p0 + p1;
      double f_ij = 4 * p0 + p1;
      af += e_ij;
      num += 255 * f_ij - e_ij * e_ij;
      ptr_mat[i] = dosage ? decode[2 * p0 + p1] : sample_from_prob_uchar(p0, p1);
    }
  }

  // https://doi.org/10.1038/nrg2796
  double coef = 255 * (2 * nona);
  *info = 1 - num * 2 * nona / (af * (coef - af));
  *freq = 1 - af / coef;

  delete[] buffer_in;
  delete[] buffer_out;

  return id;
}

/******************************************************************************/

// [[Rcpp::export]]
List read_bgen(std::string filename,
               const NumericVector& offsets,
               const Environment& BM,
               const IntegerVector& ind_row,
               const IntegerVector& ind_col,
               const RawVector& decode,
               bool dosage,
               int N,
               int ncores) {

  XPtr<FBM_RW> xpBM = BM["address_rw"];
  unsigned char* ptr_mat = static_cast<unsigned char*>(xpBM->matrix());
  std::size_t n = xpBM->nrow();

  int K = offsets.size();
  myassert_size(ind_col.size(), K);
  CharacterVector ID(K);
  std::vector<double> INFO(K, NA_REAL), FREQ(K, NA_REAL);

  #pragma omp parallel num_threads(ncores)
  {
    // connection to BGEN file
    std::ifstream stream(filename.c_str(), std::ifstream::binary);
    if (!stream) Rcpp::stop("Error while opening '%s'.", filename);

    // read variants one by one
    #pragma omp for
    for (int k = 0; k < K; k++) {
      stream.seekg(offsets[k]);
      std::size_t j = ind_col[k] - 1;
      std::string id = read_variant(&stream, ptr_mat + n * j,
                                    ind_row, decode, dosage, N,
                                    &INFO[k], &FREQ[k]);
      #pragma omp critical
      ID[k] = id;
    }

    stream.close();
  }

  return List::create(_["ID"] = ID, _["INFO"] = INFO, _["FREQ"] = FREQ);
}

/******************************************************************************/
