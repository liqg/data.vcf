#include "llib.h"
#include <Rcpp.h>

//' Resize a list of string vectors by a number vector with a same length with the list.
//' The elements of x[i] after n[i] will be discarded, 
//' however, if the length of x[i] is less than n[i], the 'fill' string will be added in the end.
//' @param x a list of string vectors.
//' @param n a number vector.
//' @param fill default value for expanded vector.
//' @return a resized list will be return.
// [[Rcpp::export]]
std::vector<std::vector<std::string> > resize_list_string (
    std::vector<std::vector<std::string> > &x,
    std::vector<int> &n,
    std::string fill = ""){
  if(x.size() != n.size()){
    Rcpp::stop("the length of n must be equal with the length of x");
  }
  
  for(int i=0; i < n.size(); i++){
    if(x[i].size() != n[i]) x[i].resize(n[i], fill);
  }
  
  return x;
}



//' Split a str by sep into ikv list.
//' @param x a string vector
//' @return An list of {i, k, v}
//'
// [[Rcpp::export]]
Rcpp::List str_to_ikv (std::vector<std::string> &x, const char sep){
  std::vector<int> iv;
  std::vector<std::string> kv;
  std::vector<std::string> vv;
  
  for(int i=0; i < x.size(); i++){
    std::vector<std::string> xi = split(x[i], sep);
    for(int j=0; j < xi.size(); j++){
      std::string k="";
      std::string v="";
      
      std::string::size_type pos = 0;
      pos = xi[j].find_first_of("=");
      if(pos == std::string::npos){
        k = xi[j];
      }else{
        k = xi[j].substr(0, pos);
        if(pos < xi[j].size()){
          v = xi[j].substr(pos+1);
        }
      }
      
      iv.emplace_back(i+1);
      kv.emplace_back(std::move(k));
      vv.emplace_back(std::move(v));
    }
  }
  return Rcpp::List::create(Rcpp::Named("i") = iv,
                            Rcpp::Named("k") = kv,
                            Rcpp::Named("v") = vv);
}



//' collapse each vector of a list by a group of indexs
//' @param x list
//' @param g list
//' @param sep string
//' @param fill string
//' @return list
//'
// [[Rcpp::export]]
std::vector<std::vector<std::string> > collapse_group(
    std::vector<std::vector<std::string> > &x,
    std::vector<std::vector<unsigned int> > &g,
    std::string sep="",
    std::string fill=""){
  std::vector<std::vector<std::string> > out(x.size());
  for(int i=0; i < x.size(); i++){ // loop x
    std::vector<std::string> out_i(g.size(), fill);
    for(int j=0; j < g.size(); j++){ // loop group
      std::string s = g[j].size() > 0 ? x[i][g[j][0]-1] : "";
      for(int k=1; k < g[j].size(); k++){
        s += sep;
        s += x[i][g[j][k]-1];
      }
      out_i[j] = s;
    }
    out[i] = out_i;
  }
  return out;
}

//' remove consecutive duplicate characters in a string
//' 
// [[Rcpp::export]]
std::vector<std::string> uniq_char(std::vector<std::string> x, std::string y=""){
  char c = y[0];
  std::vector<std::string> res(x.size());
  for(int i=0; i<x.size(); i++){
    char s[x[i].size()+1];
    int k=0;
    if(x[i].size()>1){
      s[k++] = x[i][0];
      for(int j=1; j<x[i].size(); j++){
        if(c ? (x[i][j-1] != c  || x[i][j] != c) : x[i][j] != x[i][j-1]){
          s[k++]=x[i][j];
        }
      }
      s[k] = '\0';
    }
    res[i]=std::string(s);
  }
  return res;
}
