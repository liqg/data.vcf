#include "split.h"

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::string delims = std::string(1, delim);
    tokenize(s, elems, delims);
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

std::vector<std::string> &split(const std::string &s, const std::string& delims, std::vector<std::string> &elems) {
    tokenize(s, elems, delims);
    return elems;
}

std::vector<std::string> split(const std::string &s, const std::string& delims) {
    std::vector<std::string> elems;
    return split(s, delims, elems);
}

std::vector<std::vector<std::string> > split(const std::vector<std::string>& x, const std::string& sep){
  std::vector< std::vector<std::string> > out;
  for(int i = 0; i < x.size(); i++){
    out.push_back(split(x[i], sep));
  }
  return out;
}
