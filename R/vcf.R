#' open a vcf file
#' @param path the path for vcf file
#' @return a reader object for vcf reading
open_file <- function(x) {
  stopifnot(file.exists(x))
  if(x %in% c("/dev/stdin", "")) x <- "/dev/stdin"
  else x <- normalizePath(x)
  open_file_rcpp(x)
}

#' open a vcf file (gz-compressed or plain)
#' @param file filename
#' @return an vcf object which is an environment including the file connection, meta lines, header, etc.
open_vcf <- function(file){
  if(!file.exists(file))
    stop("file not exists")
  filepath <- file
  # open file
  con <- open_file(filepath)
  # meta line
  meta <- read_vcf_header(con)
  stopifnot(length(meta) > 0)
  # header line
  header_line <- meta[length(meta)]
  meta <- meta[-length(meta)]
  if (grepl("#CHROM", header_line)) {
    header <- strsplit2(header_line, "\t", fixed=T)[[1]]
    header <- sub("^#","", header)
  } else {
    stop("can not find header")
  }

  stopifnot(length(header) >= 8)
  standard_8cols <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  if(any(header[1:8] != standard_8cols)){
    stop(paste0("header line must be like: #", paste(standard_8cols, collapse="\t")))
  }

  if(length(header) >= 9)
    stopifnot(header[9] == "FORMAT")

  vcf <- new.env(parent = emptyenv())
  vcf$filepath <- filepath
  vcf$con <- con
  vcf$meta <- meta
  vcf$header <- header
  vcf$nvar <- 0L

  class(vcf) <- c("vcf", class(vcf))

  vcf
}

split_info <- function(vars, info_keys = vector("character")) {
  info_dt <- kvsplit(vars$INFO, "=", ";")
  if (length(info_dt)>0)
    setnames(info_dt, paste0("INFO.", colnames(info_dt)))
  
  if (length(info_keys) > 0) {
    missing_cols <- setdiff(paste0("INFO.",info_keys), colnames(info_dt))
    info_dt <- cbind(info_dt, new.data.table("", nrow(info_dt), missing_cols))
  } 
  
  info_dt
}

split_format <- function(vars, format_keys = vector("character"), sample_ids=vector("character")) {
  format_dt <- NULL
  if (length(vars) > 9) {
    sample_ids <- intersect(sample_ids, colnames(vars))
    if(length(sample_ids) > 0) {
      
    }
    format_dt <- setDT(split_format_rcpp(vars$FORMAT, vars[, -(1:9), with=F]))
    if(length(format_keys) > 0) {
      tmp <- expand.grid(format_keys, colnames(vars)[10:length(vars)])
      missing_cols <- setdiff(paste("FORMAT", tmp[,1], tmp[,2], sep="."), colnames(format_dt))
      format_dt <- cbind(format_dt, new.data.table("", nrow(format_dt), missing_cols))
    }
  }
  format_dt
}

#' Read n variants from a opened vcf file.
#' @param vcf A vcf object returned by open_vcf
#' @param n The number of lines to be read
#' @param info_keys A subset of INFO keys to keep, default: empty, all will be retained.
#' @param info_keys_exclude A subset of INFO keys to discard, default: empty, none will be discarded.
#' @param format_keys A subset of FORMAT keys to keep, default: empty, all will be retained.
#' @param format_keys_exclude A subset of FORMAT keys to discard, default: empty, none will be discarded.
#' @param sample_ids A subset of sample ids to keep, default: empty, all will be retained.
#' @param break_alt If TRUE, the ALT variants will be broken, default: FALSE.
#' @param break_info_keys Which info_keys will be splited when break_alt is TRUE, default: AC and AF.
#' @param break_format_keys Which break_format_keys will be split when break_alt is TRUE, default GT and AD.
#' @param split_info If TRUE, the INFO will be splited as INFO.key, INFO.AF for example, default: FALSE.
#' @param split_format If TRUE, the FORMAT will be splited as FORMAT.key.sampleid.
#' @details read_vars returns a data.table.
read_vars <- function(
  vcf, 
  n = 0L,
  info_keys = vector("character"),
  info_keys_exclude = vector("character"),
  format_keys = vector("character"),
  format_keys_exclude = vector("character"),
  sample_ids = vector("character"), 
  break_alt = FALSE,
  break_info_keys = c("AC", "AF"), 
  break_format_keys = c("GT", "AD"),
  split_info = FALSE,
  split_format = FALSE
)
{
  info_keys <- unique(info_keys)
  format_keys <- unique(format_keys)
  break_info_keys <- unique(break_info_keys)
  break_format_keys <- unique(break_format_keys)
  
  if(is.character(vcf)){
    vcf <- open_vcf(vcf)
  }
  stopifnot(inherits(vcf, "vcf"))
  if(length(sample_ids) > 0) {
    stopifnot(length(vcf$header)>9)
    stopifnot(all(sample_ids %in% vcf$header[-(1:9)]))
    sample_ids <- which(vcf$header[-(1:9)] %in% sample_ids) - 1
  } else {
    sample_ids <- vector("integer")
  }
  cnms <- vcf$header[1:8]
  if (length(vcf$header) > 9) {
    if (length(sample_ids) == 0) cnms <- vcf$header
    else cnms <- c(cnms, "FORMAT", vcf$header[-(1:9)][sample_ids+1])
  }
  
  vars <- read_lines(vcf$con, n);
  if(length(vars) == 0) return(data.table())
  
  vars <- split_var_lines_rcpp(
    vars, info_keys, format_keys, sample_ids,
    break_alt, break_info_keys, break_format_keys, info_keys_exclude, format_keys_exclude)
  setDT(vars)
  setnames(vars, cnms)
  vars[, POS:=as.integer(POS)]
  
  info_dt <- NULL
  if (split_info) {
    info_dt <- split_info(vars, info_keys)
  }
  
  format_dt <- NULL
  if (split_format && length(vars) > 9) {
    format_dt <- split_format(vars, format_keys)
  }

  cbind(vars, info_dt, format_dt)
}

#' Convert a splited data.table to a standard VCF data.table
#'
#'
tovcf <- function(DT){
  stopifnot(inherits(DT, "data.table"))
  stopifnot(all(c("CHROM","REF","ALT","POS") %in% colnames(DT)))

  if(!"ID" %in% colnames(DT)){
    ID <- "."
  }else{
    ID <- DT$ID
  }

  if(!"QUAL" %in% colnames(DT)){
    QUAL <- "."
  }else{
    QUAL <- DT$QUAL
  }

  if(!"FILTER" %in% colnames(DT)){
    FILTER <- "."
  }else{
    FILTER <- DT$FILTER
  }

  #gather INFO
  info_cols <- grep('^INFO\\.', colnames(DT), value=T)
  if(length(info_cols)==0){
    INFO <- "."
  }else{
    dt <- DT[, info_cols, with=F]
    setnames(dt, sub("^INFO.", "", colnames(dt)))
    INFO <- kvpaste(dt, "=", ";")
  }

  res <- cbind(
    CHROM=DT$CHROM,
    POS=DT$POS,
    ID=ID,
    REF=DT$REF,
    ALT=DT$ALT,
    QUAL=QUAL,
    FILTER=FILTER,
    INFO=INFO
  )

  #gather FORMAT
  format_cols <- grep("^FORMAT\\.", colnames(DT), value=T)
  if(length(format_cols) > 0){
    format_keys <- unique(stringr::str_match(format_cols, "FORMAT.([^.]+)")[,2])
    if("GT" %in% format_keys){
      format_keys <- c("GT", format_keys[format_keys!="GT"])
    }
    FORMAT <- paste(format_keys, collapse=":")

    #gather samples
    sampleids <- unique(stringr::str_match(format_cols, "FORMAT\\.([^.]+)\\.(.+)")[,3])
    res_samples <- setDT(lapply(sampleids, function(sid){
      do.call(paste, c(lapply(format_keys, function(fk){
        col <- paste0("FORMAT.", fk, ".", sid)
        if(col %in% colnames(DT)){
          v <- DT[[col]]
          v[is.na(v)] <- "."
        }else{
          v <- "."
        }
        v
      }), sep=":"))
    }))
    setnames(res_samples, sampleids)
    res <- cbind(res, FORMAT=FORMAT, res_samples)
  }
  res
}

paste_info <- function(DT) {
  info_cols <- grep('^INFO\\.', colnames(DT), value=T)
  if(length(info_cols)==0){
    INFO <- "."
  }else{
    dt <- DT[, info_cols, with=F]
    setnames(dt, sub("^INFO.", "", colnames(dt)))
    INFO <- kvpaste(dt, "=", ";")
  }
  INFO
}
