# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' open a vcf file
#' @param path the path for vcf file
#' @return a reader object for vcf reading
open_file <- function(path) {
    .Call('_data_vcf_open_file', PACKAGE = 'data.vcf', path)
}

#' close file.
#' @param reader a pointer returned by open_file.
close_file <- function(reader) {
    invisible(.Call('_data_vcf_close_file', PACKAGE = 'data.vcf', reader))
}

#' read n lines from the file.
#' @param reader a pointer returned by open_vcf.
#' @param n how many lines to read default: n=0 for file end
#' @return a character vector, return a empty vector (length=0) when reading a file that had reached EOF; return a vector of length 1 with "" for blank line.
#' TODO: add break param to determine if break the lines
read_lines <- function(reader, n) {
    .Call('_data_vcf_read_lines', PACKAGE = 'data.vcf', reader, n)
}

#' Resize a list of string vectors by a number vector with a same length with the list.
#' The elements of x[i] after n[i] will be discarded, 
#' however, if the length of x[i] is less than n[i], the 'fill' string will be added in the end.
#' @param x a list of string vectors.
#' @n a number vector.
#' @fill default value for expanded vector.
#' @return a resized list will be return.
resize_list_string <- function(x, n, fill = "") {
    .Call('_data_vcf_resize_list_string', PACKAGE = 'data.vcf', x, n, fill)
}

#' Split a str by sep into ikv list.
#' @param x a string vector
#' @return An list of {i, k, v}
#'
str_to_ikv <- function(x, sep) {
    .Call('_data_vcf_str_to_ikv', PACKAGE = 'data.vcf', x, sep)
}

#' looply collapse each vector of a list by a group of indexs
#' @param x list
#' @param g list
#' @param sep string
#' @param fill string
#' @return list
#'
collapse_group <- function(x, g, sep = "", fill = "") {
    .Call('_data_vcf_collapse_group', PACKAGE = 'data.vcf', x, g, sep, fill)
}

#' remove consecutive duplicate characters in a string
#' 
uniq_char <- function(x, y = "") {
    .Call('_data_vcf_uniq_char', PACKAGE = 'data.vcf', x, y)
}

vcf_clean_rcpp <- function(vcf, ovcf, info_keys, format_keys) {
    invisible(.Call('_data_vcf_vcf_clean_rcpp', PACKAGE = 'data.vcf', vcf, ovcf, info_keys, format_keys))
}

