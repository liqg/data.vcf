#include "llib.h"
#include "zstr.hpp"
#include <Rcpp.h>

// return non-zero if failed
int vcf_clean(std::string vcf, std::string ovcf, std::vector<std::string> info_keys, std::vector<std::string> format_keys)
{
  int ninfo = info_keys.size();
  std::set<std::string> info_set;
  for(int i=0; i<ninfo; i++) {
    info_set.insert(info_keys[i]);
  }

  if(format_keys.size() > 0) {
    if(!in<std::string>("GT", format_keys)) {
      std::cerr << "error: GT is missing in parameter format_keys.\n";
      return 1;
    }
  }

  zstr::ifstream input(vcf);
  std::ofstream output(ovcf);
  for(std::string line; getline( input, line ); )
  {
    if(line.size()  < 1 || line[0] == '\n') continue;
    if (line[0] == '#') {
      output << line << std::endl;
      continue;
    }
    auto fields = split(line, '\t');
    int nfield = fields.size();
    if(nfield < 8) {
      std::cerr << "error: less than 8 fields.\n";
      return 1;
    }
    // INFO
    if(ninfo > 0 && fields[7].size() > 0 && fields[7] != ".") {
      std::map<std::string, std::string> info;
      std::map<std::string, std::string> info2;
      kvsplit(info, fields[7], '=', ';');
      for(auto k: info_set){
        if(info.find(k) !=  info.end()){
          info2.insert(std::make_pair(k, info[k]));
        }
      }
      fields[7].clear();
      if(info2.size() > 0)
        kvpaste(fields[7], info2, '=', ';');
      else fields[7] = ".";
      //std::cout << fields[7] << std::endl;
    }

    // FORMAT
    if(nfield > 9 && format_keys.size() > 0){
      if(fields[8].substr(0, 2) != "GT"){
        std::cerr << "error: format must start with GT.\n";
        return 1;
      }
      auto informat = split(fields[8], ':');
      std::vector<std::string> oformat;
      std::vector<int> id;
      for(int i=0; i<informat.size(); i++){
        for (auto f1: format_keys){
          if (f1 == informat[i]){
            oformat.push_back(f1);
            id.push_back(i);
            break;
          }
        }
      }

      if(oformat.size() != informat.size()){
        fields[8] = join(oformat, ':');
        for (int i=9; i< fields.size(); i++) {
          auto s = split(fields[i], ':');
          if(s.size() != informat.size()){
            std::cerr << "error: format keys are not equal.\n" << line << std::endl;
            return 1;
          }
          decltype(s) news(id.size());
          for(int i=0; i<id.size(); i++){
            news[i] = s[id[i]];
          }
          fields[i] = join(news, ':');
        }
      }
    }
    output << join(fields, '\t') << std::endl;
  }

  output.close();

  return 0;
}

// [[Rcpp::export]]
void vcf_clean_rcpp(SEXP vcf, SEXP ovcf, SEXP info_keys, SEXP format_keys)
{
  int ret = vcf_clean(Rcpp::as<std::string>(vcf),
                      Rcpp::as<std::string>(ovcf),
                      Rcpp::as<std::vector<std::string>>(info_keys),
                      Rcpp::as<std::vector<std::string>>(format_keys));
  if(ret != 0) {
    Rcpp::stop("");
  }
}
