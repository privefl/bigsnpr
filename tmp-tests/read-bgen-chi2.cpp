/******************************************************************************/

#include <bigstatsr/BMAcc.h>
#include <zlib.h>
#include "fstream.h"

using namespace Rcpp;

/******************************************************************************/

inline unsigned char sample_from_prob(double p0, double p1) {
  double first = ::unif_rand() * 255 - p0;
  return((first < 0) ? 4 : ((first < p1) ? 5 : 6));
}

/******************************************************************************/

inline double chi2(double p_obs, double p_exp, double offset) {
  double num = std::abs(p_obs - p_exp) - offset;
  return num * num / p_exp;
}

inline double hwe_chisq(double sum_pAA, double sum_pAB, double nona) {
  // https://doi.org/10.1002/gepi.20612
  double pAA = sum_pAA / 255 / nona;
  double pAB = sum_pAB / 255 / nona;
  double pBB = 1 - (pAA + pAB);
  double pA = (sum_pAA + sum_pAB / 2) / 255 / nona;
  double pB = 1 - pA;
  double offset = 0.5 / nona;
  double chi2_1 = chi2(pAA, pA * pA,     offset);
  double chi2_2 = chi2(pAB, 2 * pA * pB, offset);
  double chi2_3 = chi2(pBB, pB * pB,     offset);
  return nona * (chi2_1 + chi2_2 + chi2_3);
}

/******************************************************************************/

std::string read_variant(std::ifstream * ptr_stream,
                         unsigned char * ptr_mat,
                         const IntegerVector& ind_row,
                         const RawVector& decode,
                         bool dosage,
                         int N,
                         double * info,
                         double * freq,
                         double * hwe) {

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

  // read decompress "probabilities" and store them as rounded dosages or hard calls
  int n = ind_row.size();
  int nona = n;
  double af = 0, num = 0, sum_pAA = 0, sum_pAB = 0;
  for (int i = 0; i < n; i++) {
    int i_G = ind_row[i];
    int i_pld = 8 + i_G;
    if (buffer_out[i_pld] >= 0x80) {
      ptr_mat[i] = 3;  // missing
      nona--;
    } else {
      // probabilities * 255
      int i_prblt = 10 + N + 2 * i_G;
      unsigned char p0 = buffer_out[i_prblt];
      unsigned char p1 = buffer_out[i_prblt + 1];
      sum_pAA += p0;
      sum_pAB += p1;
      double e_ij = 2 * p0 + p1;
      double f_ij = 4 * p0 + p1;
      af += e_ij;
      num += 255 * f_ij - e_ij * e_ij;
      ptr_mat[i] = dosage ? decode[2 * p0 + p1] : sample_from_prob(p0, p1);
    }
  }

  // https://doi.org/10.1038/nrg2796
  double coef = 255 * (2 * nona);
  *info = 1 - num * 2 * nona / (af * (coef - af));
  *freq = 1 - af / coef;
  *hwe  = hwe_chisq(sum_pAA, sum_pAB, nona);

  delete[] buffer_in;
  delete[] buffer_out;

  return id;
}

/******************************************************************************/

// [[Rcpp::export]]
List read_bgen(std::string filename,
               NumericVector offsets,
               Environment BM,
               IntegerVector ind_row,
               IntegerVector ind_col,
               RawVector decode,
               bool dosage,
               int N,
               int ncores) {

  XPtr<FBM_RW> xpBM = BM["address_rw"];
  unsigned char* ptr_mat = static_cast<unsigned char*>(xpBM->matrix());
  std::size_t n = xpBM->nrow();

  int K = offsets.size();
  myassert_size(ind_col.size(), K);
  CharacterVector ID(K);
  std::vector<double> INFO(K, NA_REAL), FREQ(K, NA_REAL), HWE(K, NA_REAL);

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
                                    &INFO[k], &FREQ[k], &HWE[k]);
      #pragma omp critical
      ID[k] = id;
    }

    stream.close();
  }

  return List::create(
    _["ID"] = ID, _["INFO"] = INFO, _["FREQ"] = FREQ, _["HWE"] = HWE);
}

/******************************************************************************/
