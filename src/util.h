#ifndef _UTIL_H_79
#define _UTIL_H_79 1

#include <Rcpp.h>
#include <string>       // std::string
#include <sstream>      // std::istringstream
#include "llib.h"

static inline void getptr(void **ret, SEXP x, int i)
{
  switch (TYPEOF(x)) {
  case VECSXP:
    *ret = VECTOR_ELT(x, i);
    break;
  case INTSXP: case LGLSXP:
    *ret = INTEGER(x) + i;
    break;
  case STRSXP:
    *ret = (void*)(CHAR(STRING_ELT(x, i)));
    break;
  case REALSXP:
    *ret = REAL(x) + i;
    break;
  default:
    *ret = nullptr;
    std::cerr << "not supported type: " << Rf_type2char(TYPEOF(x)) << std::endl;
    break;
  }
}

static inline void getptr(void **ret, SEXP df, int i, int j)
{
  *ret = nullptr;
  if(TYPEOF(df) == VECSXP && Rf_length(df) > j) {
    getptr(ret, VECTOR_ELT(df, j), i);
  }
};

template<typename T>
static inline T getptr(SEXP x, int i)
{
  void *p;
  getptr(&p, x, i);
  return (T)p;
}

template<typename T>
static inline T getptr(SEXP df, int i, int j)
{
  void *p;
  getptr(&p, df, i, j);
  return (T)p;
}

template<typename T>
static inline T getptr(const Rcpp::DataFrame &df, int i, int j)
{
  void *p;
  getptr(&p, df, i, j);
  return (T)p;
}

template<typename T>
static inline T getptr(const Rcpp::DataFrame &df, int i, const std::string &jnm)
{
  void *p;
  getptr(&p, df[jnm], i);
  return (T)p;
}

template<typename V, typename T>
static inline void dfset(Rcpp::DataFrame &df, const std::string &cname, int i, T x)
{
  V col = df[cname];
  col[i] = x;
};

#endif /* _UTIL_H_79 */


