/******************************************************************************/

#include <fstream>

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
  char *buffer = new char[len + 1];
  ptr_stream->read(buffer, len);
  buffer[len] = '\0';

  std::string str(buffer, len);

  delete[] buffer;

  return str;
}

/******************************************************************************/
