#ifndef BED_ACC_H
#define BED_ACC_H

/******************************************************************************/

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/noncopyable.hpp>
#include <Rcpp.h>

using namespace boost::interprocess;
using namespace Rcpp;
using std::size_t;

/******************************************************************************/

class bedAcc {
public:
  bedAcc(const std::string path, int n, int p, const RawMatrix& decode) :
  n(n), p(p), n_byte((n + 3) / 4), _lookup_byte(decode) {

    // Memory-map the bed file
    try {
      this->file = file_mapping(path.c_str(), read_only);
    } catch(interprocess_exception& e) {
      throw std::runtime_error("File not found.");
    }
    this->file_region = mapped_region(this->file, read_only);
    this->file_data =
      static_cast<const unsigned char*>(this->file_region.get_address());

    // Verify magic numbers (http://zzz.bwh.harvard.edu/plink/binary.shtml)
    if (this->file_data[0] != '\x6C' || this->file_data[1] != '\x1B')
      stop("File is not a binary PED file (.bed).");
    if (this->file_data[2] != '\x01')
      stop("File is not in SNP-major mode.");

    // Point after this magic number
    this->file_data += 3;

    // Check if given dimensions match the file
    if ((3 + this->n_byte * this->p) != this->file_region.get_size())
      stop("The provided 'n' and 'p' do not match the dimensions of the file.");
  };

  size_t nrow() const { return n; }
  size_t ncol() const { return p; }

  inline const unsigned char operator() (size_t i, size_t j) {
    const unsigned char byte = file_data[i / 4 + j * n_byte];
    return _lookup_byte(i % 4, byte);
  }

protected:
  boost::interprocess::file_mapping file;
  boost::interprocess::mapped_region file_region;
  const unsigned char* file_data;
  size_t n;
  size_t p;
  size_t n_byte;
  RawMatrix _lookup_byte;
};

/******************************************************************************/

#endif // BED_ACC_H
