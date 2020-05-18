#' clean vcf by specified keep info_keys and format_keys.
#' @param infile input vcf file name.
#' @param outfile output vcf file name.
#' @param info_keys character vector, missing for keeping all.
#' @param format_keys character vector, missing for keeping all too.
vcf_clean <- function(infile, outfile, info_keys, format_keys){
  if(infile == "stdin") infile <- "/dev/stdin"
  if(outfile == "stdout") outfile <- "/dev/stdout"
  if(infile != "/dev/stdin") infile <- normalizePath(infile, mustWork=T)
  if(outfile != "/dev/stdout") outfile <- normalizePath(outfile, mustWork=T)
  if(missing(info_keys)) info_keys <- vector("character")
  if(missing(format_keys)) format_keys <- vector("character")
  info_keys <- na.omit(info_keys[info_keys !=""])
  format_keys <- na.omit(format_keys[format_keys !=""])
  vcf_clean_rcpp(infile, outfile, info_keys, format_keys)
}
