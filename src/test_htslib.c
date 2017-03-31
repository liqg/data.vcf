#include "htslib/hts.h"
#include "htslib/synced_bcf_reader.h"

void test_line (char *input){
  bcf_srs_t *sr =  bcf_sr_init() ;  //synced reader alloc
  bcf_sr_add_reader (sr, input ); //open the file named 'input'
  bcf1_t *line; //VCF line gets put in here
  while(bcf_sr_next_line (sr)) { //loop through file
    line =  bcf_sr_get_line(sr, 0);  //read a line
  }
}

main(int argc, char *argv[]){
  test_line(argv[0]);
  exit(0);
}