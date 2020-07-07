#include <Rcpp.h>
#include <zlib.h>
#include "zstr.hpp"

struct gz_reader {
  std::string path;
  zstr::ifstream con;
  gz_reader(std::string x) :
    path(x),
    con(x){}
};

typedef Rcpp::XPtr<gz_reader> XPtr_gz_reader;