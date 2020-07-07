#include <Rcpp.h>
#include <zlib.h>
#include <Rinternals.h> // R_ExpandFileName
#include "llib.h"
#include "data.vcf_types.h"
#include "vcf.h"
#include <Rcpp.h>

using namespace Rcpp;

//' open a vcf file
//' @param path the path for vcf file
//' @return a reader object for vcf reading
// [[Rcpp::export]]
XPtr_gz_reader open_file_rcpp(std::string path){
  gz_reader *gr =  new gz_reader(path);
  XPtr_gz_reader out(gr);
  return out;
}

//' read n lines from the file.
//' @param reader a pointer returned by open_vcf.
//' @param n how many lines to read default: n=0 for file end
//' @return a character vector, return a empty vector (length=0) when reading a file that had reached EOF; return a vector of length 1 with "" for blank line.
// [[Rcpp::export]]
std::vector<std::string> read_lines(
    XPtr_gz_reader reader,
    unsigned int n)
{
  std::vector<std::string> lines;
  int i=0;
  if (n <= 0) {
    for (std::string line; getline(reader->con, line); i++) {
      lines.push_back(line);
    }
  } else {
    lines.resize(n);
    for (; i < n && getline(reader->con, lines[i]); i++) {}
    lines.resize(i);
  }
  return(lines);
}

//[[Rcpp::export]]
std::vector<std::string> read_vcf_header(XPtr_gz_reader reader) {
  std::vector<std::string> head;
  for (std::string line; getline(reader->con, line);) {
    head.push_back(line);
    if (line.substr(0, 6) == "#CHROM") break;
  }
  return head;
}

Variant::Variant(const std::string &line)
{
  if(line.size() == 0) return;
  fields = split(line, '\t');
  
  int nfield = fields.size();
  if(nfield < 8) {
    std::cerr << "error: less than 8 fields.\n" << line << std::endl;
    return;
  }
  alt = split(fields[4], ',');
  return;
}

int Variant::split_info () {
  if (kvinfo.size() > 0 || fields.size() < 8 || fields[7] == "" || fields[7] == ".") {
    return 0;
  }
  kvsplit(kvinfo, fields[7], '=', ';');
  return 1;
}

int Variant::split_format () {
  int nfield = fields.size();
  if (format.size() == 0 && nfield > 8) {
    format = split(fields[8], ':');
    assert(var.format.size() > 0);
    if (format[0] != "GT") {
      std::cerr << "error: GT is not the first one in FORMAT.\n" << join(fields, '\t') << std::endl;
      return 0;
    }
    int nsample = nfield - 9;
    if (nsample > 0) {
      sample.resize(nsample);
      for (int i=0; i<nsample; i++) {
        if (fields[i+9] != "" && fields[i+9] != ".")
          sample[i] = split(fields[i+9], ':');
        sample[i].resize(format.size());
      }
    }
  }
  return 1;
}

void Variant::subset_info(std::set<std::string> info_keys, bool exclude) {
  if (info_keys.size() > 0) {
    split_info();
    if (kvinfo.size() > 0) {
      std::map<std::string, std::string> kvinfo2;
      for(const auto &kv: kvinfo) {
        bool flag = info_keys.find(kv.first) != info_keys.end();
        if (exclude) flag = !flag;
        if (flag) {
          kvinfo2.insert(std::make_pair(kv.first, kv.second));
        }
      }
      kvinfo = kvinfo2;
      fields[7] = kvpaste(kvinfo, '=', ';');
    }
  }
}

// sample_ids must be unique and ordered in increase
void Variant::subset_sample(std::vector<int> sample_ids) {
  if (fields.size() > 9) {
    if (sample_ids.size() > 0) {
      for(int i=0; i<sample_ids.size(); i++) {
        assert(sample_ids[i]);
        std::swap(fields[9 + i], fields[9 + sample_ids[i]]);
      }
      fields.resize(9 + sample_ids.size());
    }
  }
}

void Variant::subset_format(std::set<std::string> format_keys, bool exclude) {
  if (fields.size() < 10 || format_keys.size() == 0) return;
  split_format();

  int n1 = format.size();
  bool format_flags[n1];
  for (int i=0; i < n1; i++)
  {
    if (format_keys.find(format[i]) == format_keys.end()) {
      format_flags[i] = false;
    } else {
      format_flags[i] = true;
    }
    if (exclude) format_flags[i] = !format_flags[i];
  }

  int j = 0;
  for (int i=0; i<n1; i++)
  {
    if(format_flags[i]) {
      if (i != j ) std::swap(format[i], format[j]);
      j++;
    }
  }
  format.resize(j);
  fields[8] = join(format, ':');

  for (int k=0; k < fields.size() - 9; k++) {
    j = 0;
    for (int i=0; i < n1; i++)
    {
      if (format_flags[i]) {
        if (i != j ) std::swap(sample[k][i], sample[k][j]);
        j++;
      }
    }
    sample[k].resize(j);
    fields[k+9] = join(sample[k], ':');
  }
}

Variant Variant::pick_alt(int ialt, std::set<std::string> break_info_keys, std::set<std::string> break_format_keys) {
  assert (ialt >= 0 && ialt < alt.size());
  if (alt.size() == 1) return *this;
  Variant ivar;
  int nfield = fields.size();
  ivar.fields.resize(nfield);
  ivar.fields[0] = fields[0];
  ivar.fields[1] = fields[1];
  ivar.fields[2] = fields[2];
  ivar.fields[3] = fields[3];
  ivar.fields[4] = alt[ialt];
  ivar.fields[5] = fields[5];
  ivar.fields[6] = fields[6];

  // INFO
  split_info();
  for(const auto& kv : kvinfo) {
    std::string newv = kv.second;
    if (break_info_keys.find(kv.first) != break_info_keys.end()) {
      auto v = split(kv.second, ',');
      if (ialt >= v.size()) {
        // std::cerr << 
        //   "error: the number of alts is not matched with the key of info: "
        //   << kv.first << ". \n" <<  join(fields, '\t') 
        //   << ", this key will be skiped." << std::endl;
        continue;
      } else {
        newv = v[ialt];
      }
    }
    ivar.kvinfo.insert(std::make_pair(kv.first, newv));
  }
  if (ivar.kvinfo.size() > 0) {
    ivar.fields[7] = kvpaste(ivar.kvinfo, '=', ';');
  }

  // Genotype
  if (nfield > 9) {
    ivar.fields[8] = fields[8];
    split_format();
    int nformat = format.size();
    int nsample = nfield - 9;
    int brk_flag[nformat];
    for (int i = 0; i < nformat; i++) {
      brk_flag[i] = break_format_keys.find(format[i]) == break_format_keys.end() ? 0 : 1;
    }
    ivar.sample.resize(nsample);
    for (int i=0; i<nsample; i++) {
      //ivar.fields[i+9] = ".";
      auto gt = split(sample[i][0], '/');
      std::string delm = "/";
      if (gt.size() == 0) {
        split(sample[i][0], '|');
        if (gt.size() == 2) delm = "|";
      }
      if (gt.size() != 2) continue;
      int a1 = std::atoi(gt[0].c_str());
      int a2 = std::atoi(gt[1].c_str());
      if(a1 < 0 || a2 < 0 || a1 > alt.size() || a2 > alt.size() || (a1 != ialt+1 && a2 != ialt+1)) continue;
      ivar.sample[i].resize(nformat);
      for (int j=0; j<nformat; j++) {
        ivar.sample[i][j] = sample[i][j];
        if (brk_flag[j]) {
          if (j==0 && format[0] == "GT") {
            std::string c1 = a1 == ialt + 1 ? "1" : "0";
            std::string c2 = a2 == ialt + 1 ? "1" : "0";
            ivar.sample[i][0] = c1 + delm + c2;
          } else {
            auto tmp = split(sample[i][j], ',');
            if (tmp.size() == alt.size()+1) {
              std::string c1 = a1 == ialt + 1 ? tmp[ialt+1] : tmp[0];
              std::string c2 = a2 == ialt + 1 ? tmp[ialt+1] : tmp[0];
              ivar.sample[i][j] = c1 + "," + c2;
            } else if (tmp.size() == alt.size()) {
              ivar.sample[i][j] = tmp[ialt];
            }
          }
        }
      }
      //ivar.fields[i+9].resize(0);
      join(ivar.fields[i+9], ivar.sample[i], ':');
    }
  }
  
  return ivar;
}

std::set<std::string> vec2set(std::vector<std::string> &x) {
  std::set<std::string> s;
  for(auto k: x) {
    s.insert(k);
  }
  return s;
}

std::vector<Variant> lines_to_vars(
    std::vector<std::string> lines,
    std::set<std::string> info_keys,
    std::set<std::string> format_keys,
    std::vector<int> sample_ids,
    bool break_alt,
    std::set<std::string> break_info_keys,
    std::set<std::string> break_format_keys,
    std::set<std::string> info_exclude_keys,
    std::set<std::string> format_exclude_keys
)
{
  std::vector<Variant> vars;
  vars.reserve(lines.size());
  for(auto l: lines) {
    Variant var(l);
    if (var.fields.size() < 8) {
      Rcpp::stop("error: the number of columns < 8.\n");
    }
    var.subset_info(info_keys, false);
    var.subset_info(info_exclude_keys, true);
    var.subset_sample(sample_ids);
    var.subset_format(format_keys, false);
    var.subset_format(format_exclude_keys, true);
    
    if (break_alt && var.alt.size() > 1) {
      for (int i=0; i < var.alt.size(); i++) {
        vars.emplace_back(var.pick_alt(i, break_info_keys, break_format_keys));
      }
    } else {
      vars.emplace_back(var);
    }
  }
  return vars;
}

// [[Rcpp::export]]
Rcpp::DataFrame split_var_lines_rcpp (
    std::vector<std::string> lines,
    std::vector<std::string> info_keys,
    std::vector<std::string> format_keys,
    std::vector<int> sample_ids,
    bool break_alt,
    std::vector<std::string> break_info_keys,
    std::vector<std::string> break_format_keys,
    std::vector<std::string> info_keys_exclude,
    std::vector<std::string> format_keys_exclude
)
{
  auto info_set = vec2set(info_keys);
  auto format_set = vec2set(format_keys);
  auto break_info_set = vec2set(break_info_keys);
  auto break_format_set = vec2set(break_format_keys);
  auto format_exclude_set = vec2set(format_keys_exclude);
  auto info_exclude_set = vec2set(info_keys_exclude);

  std::vector<std::vector<std::string>> df;
  auto vars = lines_to_vars(lines, info_set, format_set, sample_ids, break_alt, break_info_set, break_format_set, info_exclude_set, format_exclude_set);
  int nr = vars.size();
  if(nr < 1) return df;
  int nc = vars[0].fields.size();
  df.resize(nc);
  for (int i=0; i<nc; i++) {
    df[i].resize(nr);
    for (int j=0; j<nr; j++)
    {
      df[i][j] = vars[j].fields[i];
    }
  }
  
  auto ret = Rcpp::DataFrame::create(df, Rcpp::Named("stringsAsFactors") = false);

  return ret;
}

//[[Rcpp::export]]
Rcpp::DataFrame split_format_rcpp(
    Rcpp::CharacterVector format,
    Rcpp::DataFrame samples
    )
{
  std::vector<std::vector<std::string>> df;
  std::map<std::string, int> format_map;
  std::vector<std::vector<int>> format_ids(format.size());
  std::vector<std::string> format_nms;
  int format_count = 0;
  for (int i=0; i<format.size(); i++) {
    auto keys = split(std::string(format[i]), ':');
    format_ids[i].resize(keys.size());
    for (int j=0; j<keys.size(); j++) {
      auto search = format_map.find(keys[j]);
      if (search == format_map.end()) {
        format_map.emplace(std::make_pair(keys[j], format_count));
        format_nms.push_back(keys[j]);
        format_ids[i][j] = format_count;
        format_count++;
      } else {
        format_ids[i][j] = search->second;
      }
    }
  }
  df.resize(format_count*samples.size());
  int nr = samples.nrow();
  for (int i =0; i < df.size(); i++) {
    df[i].resize(nr);
  }

  for (int i=0; i<samples.size(); i++) {
    std::vector<std::string> sp = Rcpp::as<std::vector<std::string>>(samples[i]);
    for (int j=0; j<nr; j++) {
      auto s = split(sp[j], ':');
      if(s.size() != format_ids[j].size()) continue;
      for (int k=0; k<s.size(); k++) {
        df[i*format_count + format_ids[j][k]][j] = s[k];
      }
    }
  }
  
  Rcpp::DataFrame ret = Rcpp::DataFrame::create(df, Rcpp::Named("stringsAsFactors") = false);
  Rcpp::CharacterVector cnms(samples.size() * format_count);
  auto sample_nms = Rcpp::as<std::vector<std::string>>(samples.names());
  for (int i=0; i<samples.size(); i++) {
    for (int j=0; j<format_count; j++){
      cnms[i*format_count + j] ="FORMAT." + format_nms[j] + "." + sample_nms[i];
    }
  }
  ret.names() = cnms;
  
  return ret;
}