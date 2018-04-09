#include <Rcpp.h>
#include <zlib.h>

/* Size of the block of memory to use for reading. */
#define LENGTH 0x1000 // hexadecimel for 4096.

struct gz_reader {
  const char *path;
  gzFile con;
  std::string buffer; // the remaining 
};

inline void finalizer_reader(gz_reader * reader) {
  if (reader){
    gzclose(reader->con);
    free(reader);
  }
}

typedef Rcpp::XPtr<gz_reader, Rcpp::PreserveStorage, finalizer_reader> XPtr_gz_reader;