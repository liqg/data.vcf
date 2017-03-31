#include <Rcpp.h>
#include <zlib.h>
#include "split.h"
#include <Rinternals.h> // R_ExpandFileName
#include <string.h>

using namespace Rcpp;

/* Size of the block of memory to use for reading. */
#define LENGTH 0x1000 // hexadecimel for 4096.

struct ikv{
  IntegerVector i; //varid
  CharacterVector k; // key
  CharacterVector v; // value
};

struct gz_reader {
  const char *path;
  gzFile con;
  std::string remains; // the remaining 
};

inline void finalizer_reader(gz_reader * reader) {
  if (reader){
    gzclose(reader->con);
    free(reader);
  }
}
typedef Rcpp::XPtr<gz_reader, Rcpp::PreserveStorage, finalizer_reader> XPtr_gz_reader;

//' open a vcf file
//' @param path the path for vcf file 
//' @return a reader object for vcf reading
//' @export
// [[Rcpp::export]]
XPtr_gz_reader open_file(std::string path){
  gz_reader *gr =  new gz_reader;
  gr->remains = "";
  
  gr->path = R_ExpandFileName(path.c_str());
  
  gzFile file = gzopen(gr->path, "r");
  
  if(! file){
    stop("gzopen failed");
  }else{
    gr->con = file;
    gzbuffer(file, 64*1024);
  }

  XPtr_gz_reader out(gr);
  
  return out;
}

void close_file(XPtr_gz_reader reader){
  gzclose(reader->con);
}

//' read n lines from the file.
//' @param reader a pointer returned by open_vcf.
//' @param n how many lines to read default: n=0 for file end
//' @value a character vector, return a vector of length 0 when reading a file that had reached EOF; return a vector of length 1 with "" for blank line.
//' @export a character vector of <= n lines.
// [[Rcpp::export]]
std::vector<std::string> read_lines( XPtr_gz_reader reader, unsigned int n){
  std::vector<std::string> lines;
  
  while(lines.size() < n ){
    if(gzeof(reader->con) & reader->remains.size() == 0){
      break;
    }
    
    Rcpp::checkUserInterrupt();
    
    int bytes_read;
    char buffer[LENGTH];
    bytes_read = gzread (reader->con, buffer, LENGTH - 1);
    buffer[bytes_read] = '\0';
    
    if(bytes_read < LENGTH -1 & !gzeof (reader->con)){
      // stop if gzerror
      int err;
      const char * error_string;
      error_string = gzerror (reader->con, & err);
      if (err) {
        Rcpp::Rcerr << "Error: " << error_string << ".\n";
        stop("gzread err");
      }
    }
    
    // paste remains and buffer
    std::string remainder_buffer = reader->remains + std::string(buffer);

    // split lines
    char * p1 = &remainder_buffer[0]; // forward ptr
    char * p2 = &remainder_buffer[0]; // ptr for the str start
    
    while(*p1 != '\0' & lines.size() < n){
      if(*p1 == '\n'){
        *p1 = '\0';
        lines.push_back(std::string(p2));
        p2 = p1;
        p2++;
      }
      p1++;
    }
    
    reader->remains = std::string(p2);
  }
  return(lines);
}


