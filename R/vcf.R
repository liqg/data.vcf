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

#' Split sample col by FORMAT in a split-resize strategy way which is very fast.
#' The elements of x[i] after n[i] will be discarded, 
#' however, if the length of x[i] is less than n[i], the '' string will be added in the end.
#' @param x the char vector of a sample col.
#' @param FORMAT_length a length vector of FORMAT feilds.
#' @return a char vector.
reformat_sample <- function(x, FORMAT_length){
  xs <- strsplit(x,":",fixed=T)
  # ensure that the length of splited sample col is same as the splited format col.
  if(!all(sapply(xs, length) == FORMAT_length)){
    xs <- resize_list_string(xs, FORMAT_length, fill="")
  }
  unlist(xs)
}

#' reformat INFO and FORMAT as varid, key, value tables
#' @param lines data.table returned by read_body
#' @return list a list of data.tables, including C7: cols1-7, INFO:, FORMAT: FORMAT and samples' genotype information.
reformat_body <- function(lines, varid_offset=0){
  variants <- list()

  if(nrow(lines) == 0)
    return(variants)

  variants$C7 <- lines[, .(varid=.I, CHROM, POS, ID, REF, ALT, QUAL, FILTER)]

  variants$INFO <- setDT(str_to_ikv(lines$INFO, sep=";"))
  setnames(variants$INFO, "i", "varid")
  variants$INFO <- variants$INFO[!(varid %in% lines[, .I[INFO %in% c(".", "", NA)]])] # remove missing info
  
  # split sample format
  if(ncol(lines) >= 9){
    splited_format <- strsplit(lines$FORMAT, ":", fixed=TRUE)
    FORMAT_length <- sapply(splited_format, length)
    if(any(FORMAT_length == 0)){
      stop("Empty string is not allowed in the FORMAT col.")
    }
    variants$FORMAT <- data.table(
      varid = rep(1:nrow(lines), FORMAT_length),
      k = unlist(splited_format))
  }

  if(ncol(lines) >= 10){
    variants$FORMAT <- cbind(
      variants$FORMAT, 
      lines[, lapply(.SD, reformat_sample, FORMAT_length = ..("FORMAT_length")), .SDcols=10:ncol(..("lines"))])
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
table_variants <- function(variants, INFO=NULL, FORMAT=NULL){
  res <- data.table()
  if(length(variants)==0) return(res)
  ikv_info <- variants$INFO
  if(!is.null(INFO)){
    ikv_info <- variants$INFO[k %in% INFO]
  }
  info_dt <- dcast(variants$INFO, varid ~ k, value.var = "v")
  setnames(info_dt, c("varid", paste0("INFO.", setdiff(colnames(info_dt),"varid"))))
  setkey(info_dt, "varid")
  dt <- merge(variants$C7, info_dt, by="varid", all.x=T, suffixes=c("", "."))
  if("FORMAT" %in% names(variants)){
    format_ikv <- variants$FORMAT
    if(!is.null(FORMAT)){
      format_ikv <-  format_ikv[k %in% FORMAT]
    }
    unique(format_ikv, by=c("varid","k"))
    format_dt <- dcast(format_ikv, varid ~ k, value.var = setdiff(colnames(variants$FORMAT), c("varid","k")))
    newnames <- setdiff(colnames(format_dt), "varid")
    
    setnames(gt_dt, c())
  }
  
  
}

table_vcf <-function(infile, outfile, fields, n=1000L){
  stopifnot(length(fields)>0)
  # stopifnot(length(fields) == length(fieldClasses))
  # stopifnot(all(fieldClasses %in% c("numeric","integer","character")))
  newnames <- fields
  if(!is.null(names(fields))){
    newnames[names(fields) != ""] <- names(fields)[names(fields)!=""]
  }
  
  dt0 <- fread(paste0(paste(rep(1, length(fields)), collapse="\t"), "\n"),
               header=F, sep="\t", col.names=fields)[0]
  
  vcf <- open_vcf(infile);
  write(paste0(paste(rep(1, length(newnames)), collapse="\t"), "\n"), file=outfile)
  while(nrow(lines_dt <- read_body(vcf, n=n)) > 0){
    if(nrow(lines_dt)==0){
      return(NULL)
    }
    variants <- reformat_body(lines_dt)
    info_dt <- dcast(variants$INFO, varid ~ k, value.var = "v")
    setkey(info_dt, "varid")
    dt <- merge(variants$C7, info_dt, by="varid", all.x=T, suffixes=c("", "."))
    dt <- dt[, fields, with=FALSE]
    dt <- rbind(dt0, dt, fill=TRUE)
    fwrite(dt, outfile, append=T, sep="\t", col.names=FALSE)
  }
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


