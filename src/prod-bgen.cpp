/******************************************************************************/

#include <bigstatsr/arma-strict-R-headers.h>
#include <bigstatsr/BMAcc.h>
#include <zlib.h>
#include "fstream.h"

using namespace Rcpp;

/******************************************************************************/

inline double sample_from_prob(double p0, double p1) {
  double first = ::unif_rand() * 255 - p0;
  return((first < 0) ? 0 : ((first < p1) ? 1 : 2));
}

/******************************************************************************/

void read_variant(std::ifstream * ptr_stream,
                  arma::mat& X, int j,
                  const IntegerVector& ind_row,
                  const NumericVector& decode,
                  bool dosage,
                  int N) {

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

  // read uncompressed "probabilities" and compute products
  int n = ind_row.size();
  for (int i = 0; i < n; i++) {
    int i_G = ind_row[i];
    int i_pld = 8 + i_G;
    if (buffer_out[i_pld] >= 0x80) {
      X(i, j) = NA_REAL;  // missing
    } else {
      int i_prblt = 10 + N + 2 * i_G;
      // probabilities * 255
      unsigned char p0 = buffer_out[i_prblt];
      unsigned char p1 = buffer_out[i_prblt + 1];
      X(i, j) = dosage ? decode[2 * p0 + p1] : sample_from_prob(p0, p1);
    }
  }

  delete[] buffer_in;
  delete[] buffer_out;
}

/******************************************************************************/

// [[Rcpp::export]]
arma::mat& extract_submat_bgen(std::string filename,
                               const std::vector<size_t>& offsets,
                               arma::mat& X,
                               const IntegerVector& ind_row,
                               const NumericVector& decode,
                               bool dosage,
                               int N,
                               int ncores) {

  int m = offsets.size();

  #pragma omp parallel num_threads(ncores)
  {
    // connection to BGEN file
    std::ifstream stream(filename.c_str(), std::ifstream::binary);
    if (!stream) Rcpp::stop("Error while opening '%s'.", filename);

    // read variants into X
    #pragma omp for
    for (int j = 0; j < m; j++) {
      stream.seekg(offsets[j]);
      read_variant(&stream, X, j, ind_row, decode, dosage, N);
    }

    stream.close();
  }

  return X;
}

/******************************************************************************/

// [[Rcpp::export]]
arma::mat& prod_bgen2(std::string filename,
                      const NumericVector& offsets,
                      arma::mat& XY,
                      const arma::mat& Y,
                      const IntegerVector& ind_row,
                      const NumericVector& decode,
                      bool dosage,
                      int N,
                      int max_size,
                      int ncores) {

  int n = ind_row.size();
  int m = offsets.size();

  arma::mat X(n, max_size);  // temporary matrix to access blocks
  std::vector<size_t> sub_offsets(max_size);

  for (int j = 0; j < m; ) {

    int k;
    for (k = 0; k < max_size && j < m; k++, j++) sub_offsets[k] = offsets[j];

    if (k < max_size) {
      // last block can be shorter
      sub_offsets.resize(k);
      XY += extract_submat_bgen(filename, sub_offsets, X, ind_row,
                                decode, dosage, N, ncores).head_cols(k) *
        Y.rows(j - k, j - 1);
    } else {
      // k == max_size
      XY += extract_submat_bgen(filename, sub_offsets, X, ind_row,
                                decode, dosage, N, ncores) *
        Y.rows(j - k, j - 1);
    }
  }

  return XY;
}

/******************************************************************************/
