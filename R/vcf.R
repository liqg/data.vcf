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
  meta <- c()
  while(length(l <- read_lines(con, 1)) > 0 & substr(l, 1, 2) == "##"){
    meta <- c(meta, l)
  }
  # header line
  if(substr(l, 1, 4) == "#CHR"){
    header <- strsplit(l, "\t", fixed=T)[[1]]
    header <- sub("^#","", header)
  }else{
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

  class(vcf) <- c("vcf", class(vcf))

  vcf
}

#' read n variants
#' @param vcf a vcf object returned by open_vcf.
#' @param n n variants/lines to read
#' @return data.table the vcf table format.
read_body <- function(vcf, n=0){
  lines <- read_lines(vcf$con, n)

  if(length(lines) == 0){
    return(data.table())
  }

  lines <- tstrsplit(lines, split="\t", fixed=TRUE)
  if(length(vcf$header) != length(lines)){
    # if(length(vcf$header) == 9 & length(lines) == 8){
    #   vcf$header <- vcf$header[1:8] # VEP will get extra FORMAT
    # }else{
    #   stop("The number of cols in the body must be equal to that in the header line!")
    # }
    stop("The number of cols in the body must be equal to that in the header line!")
  }
  names(lines) <- vcf$header
  setDT(lines)
  lines
}

#' reformat INFO and FORMAT as varid, key, value tables
#' @param lines data.table returned by read_body
#' @return list a list of data.tables, including C7: cols1-7, INFO:, FORMAT: FORMAT and samples' genotype information.
reformat_body <- function(lines, varid_offset=0){
  variants <- list()

  if(nrow(lines) == 0)
    return(variants)

  variants$C7 <- lines[, .(varid=.I, CHROM, POS, ID, REF, ALT, QUAL, FILTER)]

  variants$INFO <- setDT(split_info(lines$INFO))
  setnames(variants$INFO, "i", "varid")
  variants$INFO <- variants$INFO[!(varid %in% lines[, .I[INFO %in% c(".", "", NA)]])] # remove missing info

  if(ncol(lines) >= 9){
    ans <- strsplit(lines$FORMAT, ":", fixed=TRUE)
    FORMAT_Length <- sapply(ans, length)
    if(any(FORMAT_Length == 0)){
      stop("Empty string is not allowed in the FORMAT col.")
    }
    variants$FORMAT <- data.table(
      varid = rep(1:nrow(lines), FORMAT_Length),
      k = unlist(ans))
  }

  if(ncol(lines) >= 10){
    samples_split <- lines[, 10:ncol(lines)][
      ,
      lapply(.SD, function(x){
        gt <- strsplit(x,":",fixed=T)
        if(!all(sapply(gt, length) == FORMAT_Length)){
          gt <- resize_list_string(gt, FORMAT_Length)
        }
        unlist(gt)
      })]
    variants$FORMAT <- cbind(variants$FORMAT, samples_split)
  }

  variants$C7[, varid:=varid+varid_offset]
  setkey(variants$C7, "varid")
  variants$INFO[, varid:=varid+varid_offset]
  setkey(variants$INFO, "varid")
  if("FORMAT" %in% names(variants)){
    variants$FORMAT[, varid:=varid+varid_offset]
    setkey(variants$FORMAT, "varid")
  }
  variants
}

#' revert variants back to the vcf body table format.
#' @param variants an object returned by reformat_body.
#' @return data.table
tobody <- function(variants){
  if(!all(c("C7","INFO") %in% names(variants))){
    stop("wrong variant object")
  }else{
    setkey(variants$C7, "varid")
    setkey(variants$INFO, "varid")
  }

  if("FORMAT" %in% names(variants)){
    setkey(variants$FORMAT, "varid")
  }

  lines <- copy(variants$C7)

  variants$INFO[, kv:=paste(k, v, sep="=")]
  variants$INFO[is.na(v), kv:=paste(k, ".", sep="=")] # reset missing values
  variants$INFO[v == "", kv:=as.character(k)] # reset flag info
  info_dt <- variants$INFO[, .(INFO=paste(kv, collapse=";")), by="varid"]
  lines[info_dt, INFO:=INFO]
  lines[is.na(INFO), INFO:="."]
  variants$INFO[, kv:=NULL]

  if("FORMAT" %in% names(variants)){
    if(nrow(variants$FORMAT) > 0){
      inds <- setdiff(colnames(variants$FORMAT), c("k","varid"))
      varid_groups <- tapply(1:nrow(variants$FORMAT), variants$FORMAT$varid, c)
      format_dt <- collapse_group(variants$FORMAT[, c("k", inds), with=FALSE], varid_groups, sep=":")
      names(format_dt) <- c("FORMAT", inds)
      format_dt$varid <- as.integer(names(varid_groups))
      setDT(format_dt, key="varid")
      lines <- merge(lines, format_dt[, c("varid","FORMAT", inds), with=F], all.x=TRUE)
    }
  }
  lines[, -"varid"]
}

