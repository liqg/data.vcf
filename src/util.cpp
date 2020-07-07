#include <Rcpp.h>
#include "util.h"

// [[Rcpp::export]]
std::vector<std::string> kvpaste_rcpp(Rcpp::DataFrame input, std::string sep1, std::string sep2, bool nakey){
  std::vector<std::string> out(input.nrow());
  auto keys = Rcpp::as<std::vector<std::string>>(input.names());
  if (input.size() < 1) return out;

  for (int j=0; j<input.ncol(); j++) {
    Rcpp::CharacterVector col = Rcpp::as<Rcpp::CharacterVector>(input[j]);
    for (int i=0; i<input.nrow(); i++) {
      auto v = col[i];
      if(v != "") {
        if(v != NA_STRING) {
          if(out[i].size() > 0) out[i] += sep2;
          out[i] = out[i] + keys[j] + sep1 + v ;
        }
        else if (nakey) {
          if(out[i].size() > 0) out[i] += sep2;
          out[i] = out[i] + keys[j];
        }
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
int count_beyond_threshold(Rcpp::NumericVector x, int y) {
  // count how many non-dup elements in vector x  that are >= y times
  // x=1,2,2,3,3,3,4; y=3
  // 1
  std::sort(x.begin(), x.end());
  int xn = x.size();
  if(y > xn){
    return 0;
  }

  int element = x[0];
  int nelement = 1;
  int res = 0;
  if(nelement == y){
    res++;
  }

  for(int i=1; i < xn; i++){
    if(x[i] == element){
      nelement++;
    }else{
      element = x[i];
      nelement = 1;
    }
    if(nelement == y){
      res++;
    }
  }
  return res;
}

//[[Rcpp::export]]
Rcpp::DataFrame kvsplit_rcpp(std::vector<std::string> input, char sep1, char sep2, bool na){
  std::vector<Rcpp::CharacterVector> out;
  Rcpp::CharacterVector cnms;
  std::map<std::string, int> idmap;
  int nr = input.size();
  for(int i=0; i<nr; i++){
    std::map<std::string, std::string> attr;
    kvsplit(attr, input[i], sep1, sep2);
    for (const auto &kv: attr) {
      auto id = idmap.find(kv.first);
      if (idmap.find(kv.first) == idmap.end()) {
        idmap.insert(std::make_pair(kv.first, idmap.size()));
        Rcpp::CharacterVector vec(nr);
        if(na) vec.fill(NA_STRING);
        vec[i] = kv.second;
        out.emplace_back(vec);
        cnms.push_back(kv.first);
      } else {
        out[id->second][i] = kv.second;
      }
    }
  }
  Rcpp::DataFrame ret = Rcpp::DataFrame::create(out, Rcpp::Named("stringsAsFactors") = false);
  ret.names() = cnms;
  return ret;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix col_sums(
    Rcpp::NumericMatrix m,
    Rcpp::List rows_list
){
  int ncols = m.ncol();
  int nlist = rows_list.size();
  Rcpp::NumericMatrix res(nlist, ncols);
  
  for (int i = 0; i < nlist; i++){
    Rcpp::NumericVector rows = Rcpp::as<Rcpp::NumericVector>(rows_list[i]) - 1; //0-based
    int nr = rows.size();
    if(nr == 0){
      for ( int j = 0; j < ncols; j++){
        res(i, j) = 0;
      }
    }
    else{
      for ( int j = 0; j < ncols; j++){
        double sum_j = m(rows[0], j);
        if(nr > 1){
          for(int k = 1; k < nr; k++){
            sum_j = sum_j + m(rows[k], j);
          }
        }
        res(i, j) = sum_j;
      }
    }
  }
  return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix col_maxs(
    Rcpp::NumericMatrix m,
    Rcpp::List rows_list
){
  int ncols = m.ncol();
  int nlist = rows_list.size();
  Rcpp::NumericMatrix res(nlist, ncols);

  for (int i = 0; i < nlist; i++){
    Rcpp::NumericVector rows = Rcpp::as<Rcpp::NumericVector>(rows_list[i]) - 1; //0-based
    int nr = rows.size();
    if(nr == 0){
      for ( int j = 0; j < ncols; j++){
        res(i, j) = 0;
      }
    }
    else{
      for ( int j = 0; j < ncols; j++){
        double max_j = m(rows[0], j);
        if(nr > 1){
          for(int k = 1; k < nr; k++){
            max_j = max_j > m(rows[k], j)? max_j : m(rows[k], j);
            // max_j = max_j + m(rows[k], j);
          }
        }
        res(i, j) = max_j;
      }
    }
  }
  return Rcpp::wrap(res);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// col_maxs(m=matrix(1:12,nrow=4), list(1,2,1:2))

