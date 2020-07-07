#ifndef VCF_H_7_9
#define VCF_H_7_9
#include "llib.h"

class Variant {
public:
  std::vector<std::string> fields;
  std::vector<std::string> alt;
  std::map<std::string, std::string> kvinfo;
  std::vector<std::string> format;
  std::vector<std::vector<std::string>> sample;
  Variant(){};
  Variant(const std::string &);
  int split_info ();
  int split_format ();
  Variant pick_alt(int ialt, std::set<std::string> break_info_keys, std::set<std::string> break_format_keys);
  void subset_info(std::set<std::string> info_keys, bool exclude);
  void subset_sample(std::vector<int> sample_ids); // sample_ids must be ordered in increase
  void subset_format(std::set<std::string> format_keys, bool exclude);
};


#endif