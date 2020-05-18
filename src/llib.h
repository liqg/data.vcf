#ifndef LLIB_H_79
#define LLIB_H_79 1

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <set>
#include <map>

// from Marius, http://stackoverflow.com/a/1493195/238609
// thanks to Evan Teran, http://stackoverflow.com/questions/236129/how-to-split-a-string/236803#236803

template < class ContainerT >
static void tokenize(const std::string& str, ContainerT& tokens,
              const std::string& delimiters = " ", const bool trimEmpty = false)
{
  std::string::size_type pos, lastPos = 0;
  while(true)
  {
    pos = str.find_first_of(delimiters, lastPos);
    if(pos == std::string::npos)
    {
      pos = str.length();
      if(pos != lastPos || !trimEmpty) {
        tokens.push_back(typename ContainerT::value_type(str.data()+lastPos, (typename ContainerT::value_type::size_type)pos-lastPos));
      }
      break;
    }
    else
    {
      if(pos != lastPos || !trimEmpty) {
        tokens.push_back(typename ContainerT::value_type(str.data()+lastPos, (typename ContainerT::value_type::size_type)pos-lastPos));
      }
    }
    lastPos = pos + 1;
  }
};

static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::string delims = std::string(1, delim);
  tokenize(s, elems, delims);
  return elems;
}

static std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  return split(s, delim, elems);
}

static std::vector<std::string> &split(const std::string &s, const std::string& delims, std::vector<std::string> &elems) {
  tokenize(s, elems, delims);
  return elems;
}

static std::vector<std::string> split(const std::string &s, const std::string& delims) {
  std::vector<std::string> elems;
  return split(s, delims, elems);
}

static std::vector<std::vector<std::string>> split(const std::vector<std::string>& x, const std::string& sep){
  std::vector< std::vector<std::string>> out(x.size());
  for(int i = 0; i < x.size(); i++){
    out[i]=split(x[i], sep);
  }
  return out;
}

template <typename T>
static inline std::set<T> vec2set(std::vector<T> &x){
  std::set<T> out;
  for(auto k: x) {
    if(out.find(k) != out.end()) out.insert(k);
  }
  return out;
}

template <typename T>
static inline bool in(T x, std::vector<T> &y){
  for(auto k: y){
    if(k == x) return true;
  }
  return false;
}

template <typename T>
static inline bool in(T x, std::set<T> &y){
  return y.find(x) != y.end();
}

template <typename T>
static std::vector<T> unique(std::vector<T> &x)
{
  std::vector<T> out;
  std::set<T> xmap;
  for(auto k : x){
    if(xmap.find(k) == xmap.end()){
      xmap.insert(k);
      out.push_back(k);
    }
  }
  return out;
}

static void kvpaste(std::string &out, std::map<std::string, std::string> &x, char sep1, char sep2)
{
  for(auto kv: x){
    if(out.size() != 0)  out += sep2;
    out += kv.first;
    if(kv.second.size() != 0) out = out + sep1 + kv.second;
  }
  //std::cout << out << std::endl;
}

static std::string kvpaste(std::map<std::string, std::string> &x, char sep1, char sep2)
{
  std::string out;
  kvpaste(out, x, sep1, sep2);
  return out;
}

static void join(std::string &out, std::vector<std::string> &x, char sep)
{
  for(auto k: x){
    if(out.size() != 0) out += sep;
    out += k;
  }
}

static std::string join(std::vector<std::string> &x, char sep)
{
  std::string out;
  join(out, x, sep);
  return out;
}

static void kvsplit(std::map<std::string, std::string> &attr, std::string &x, char sep1, char sep2)
{
  std::vector<std::string> v = split(x, sep2);
  for (int j=0; j< (int)v.size(); j++) {
    std::string kv = v[j];
    int kstart=-1, kend=-1, vstart=-1, vend=-1, pos=-1;
    while(++pos < (int) kv.size()) {
      if(kstart == -1) {
        if(kv[pos] == ' ') continue; //skip empty before k
        else  kstart = pos;
      }

      if(kend == -1 && kstart != -1) {
        if (kv[pos] == sep1) {
          kend = pos - 1;
        }
        continue;
      }

      if(kend != -1 && vstart == -1) {
        if(kv[pos] == ' ') continue; // skip empty after sep1
        else {
          vstart = pos;
          break;
        }
      }
    }
    //std::cout << kstart << " " << kend << " " << vstart << " " << vend << std::endl;

    if (kstart != -1 && kend == -1) { // not found sep1
      pos = (int) kv.size();
      while(--pos >= kstart) {
        if (kv[pos] == ' ') continue;
        else {kend = pos; break;}
      }
    } else if (vstart != -1) {
      pos = (int) kv.size();
      while(--pos >= vstart) {
        if (kv[pos] == ' ') continue;
        else {vend = pos; break;}
      }
    }
    //std::cout << kstart << " " << kend << " " << vstart << " " << vend << std::endl;
    std::string k, v;
    if (kend == -1 ) continue;
    else k = kv.substr(kstart, kend-kstart+1);
    if (vend != -1) v = kv.substr(vstart, vend - vstart + 1);
    k.shrink_to_fit();
    v.shrink_to_fit();
    attr.insert(std::make_pair(k, v));
  }
}

static std::map<std::string, std::string> kvsplit(std::string &x, char sep1, char sep2)
{
  std::map<std::string, std::string> out;
  kvsplit(out, x, sep1, sep2);
  return out;
}

#endif

