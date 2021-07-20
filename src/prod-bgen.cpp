/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/BMAcc.h>
#include <zlib.h>
#include "fstream.h"

using namespace Rcpp;


#if defined(_OPENMP)
#include <omp.h>
#else
inline int omp_get_thread_num() { return 0; }
#endif

/******************************************************************************/

inline double sample_from_prob(double p0, double p1) {
  double first = ::unif_rand() * 255 - p0;
  return((first < 0) ? 0 : ((first < p1) ? 1 : 2));
}

/******************************************************************************/

void prod_variant(std::ifstream * ptr_stream,
                  const arma::vec& beta_trans_j,
                  arma::mat& res,
                  const IntegerVector& ind_row,
                  const NumericVector& decode,
                  bool dosage,
                  int N,
                  int L) {

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

  // read uncompressed "probabilities" and compute products
  int n = ind_row.size();
  for (int i = 0; i < n; i++) {
    int i_G = ind_row[i];
    int i_pld = 8 + i_G;
    if (buffer_out[i_pld] >= 0x80) {
      res.row(i).fill(NA_REAL);  // missing
    } else {
      // probabilities * 255
      int i_prblt = 10 + N + 2 * i_G;
      unsigned char p0 = buffer_out[i_prblt];
      unsigned char p1 = buffer_out[i_prblt + 1];
      double g_ij = dosage ? decode[2 * p0 + p1] : sample_from_prob(p0, p1);
      for (int l = 0; l < L; l++) res(i, l) += g_ij * beta_trans_j[l];
    }
  }

  delete[] buffer_in;
  delete[] buffer_out;
}

/******************************************************************************/

// [[Rcpp::export]]
arma::cube prod_bgen(std::string filename,
                     const NumericVector& offsets,
                     const arma::mat& beta_trans,
                     const IntegerVector& ind_row,
                     const NumericVector& decode,
                     bool dosage,
                     int N,
                     int ncores) {

  int m = offsets.size();
  myassert_size(beta_trans.n_cols, m);

  int L = beta_trans.n_rows;
  arma::cube res(ind_row.size(), L, ncores, arma::fill::zeros);

  #pragma omp parallel num_threads(ncores)
  {
    // connection to BGEN file
    std::ifstream stream(filename.c_str(), std::ifstream::binary);
    if (!stream) Rcpp::stop("Error while opening '%s'.", filename);

    int id = omp_get_thread_num();

    // compute (increment) product variant by variant
    #pragma omp for
    for (int j = 0; j < m; j++) {
      stream.seekg(offsets[j]);
      prod_variant(&stream, beta_trans.col(j), res.slice(id),
                   ind_row, decode, dosage, N, L);
    }

    stream.close();
  }

  return res;
}

/******************************************************************************/
