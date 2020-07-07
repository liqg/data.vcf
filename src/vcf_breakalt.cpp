#include <Rcpp.h>
#include "vcf.h"
#include "zstr.hpp"

int vcf_breakalt(
    std::string vcf, 
    std::string ovcf, 
    std::vector<std::string> break_info_keys, 
    std::vector<std::string> break_format_keys)
{
  std::set<std::string> info_set;
  for(int i=0; i<break_info_keys.size(); i++) {
    info_set.insert(break_info_keys[i]);
  }
  
  if(break_format_keys.size() > 0) {
    if(!in<std::string>("GT", break_format_keys)) {
      std::cerr << "error: GT is missing in parameter format_keys.\n";
      return 1;
    }
  }
  std::set<std::string> format_set;
  for(int i=0; i<break_format_keys.size(); i++) {
    format_set.insert(break_format_keys[i]);
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
    Variant var(line);
    int nfield = var.fields.size();
    if(nfield < 8) {
      std::cerr << "error: less than 8 fields.\n";
      return 1;
    }
    
    int nalt = var.alt.size();
    if (nalt <= 1) {
      output << line << std::endl;
    } else {
      for (int ialt=0; ialt < nalt; ialt++) {
        Variant ivar = var.pick_alt(ialt, info_set, format_set);
        output << join(ivar.fields, '\t') << std::endl;
      }
    }
  }
  return 0;
}

// [[Rcpp::export]]
void vcf_breakalt_rcpp(SEXP vcf, SEXP ovcf, SEXP info_keys, SEXP format_keys)
{
  int ret = vcf_breakalt(Rcpp::as<std::string>(vcf),
                         Rcpp::as<std::string>(ovcf),
                         Rcpp::as<std::vector<std::string>>(info_keys),
                         Rcpp::as<std::vector<std::string>>(format_keys));
  if(ret != 0) {
    Rcpp::stop("error");
  }
}
