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
  while(length(l <- read_lines(con, 1000)) > 0 & all(substr(l, 1, 2) == "##")){
    meta <- c(meta, l)
  }
  
  # header line
  headid = which(substr(l, 1, 4) == "#CHR")[1]
  if (!is.na(headid)) {
    header <- strsplit2(l[headid], "\t", fixed=T)[[1]]
    header <- sub("^#","", header)
  }else{
    stop("can not find header")
  }
  if(headid > 1) {
    meta <- c(meta, l[1:(headid-1)])
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
  vcf$left <- vector(mode="character")
  if (headid < length(l)) vcf$left <- l[(headid+1):length(l)]

  class(vcf) <- c("vcf", class(vcf))

  vcf
}

close.vcf <- function(vcf){
  close_file(vcf$con)
}

#' read n variants
#' @param vcf a vcf object returned by open_vcf.
#' @param n n variants/lines to read
#' @return data.table the vcf table format.
read_body <- function(vcf, n=0){
  if(length(vcf$left) > 0) {
    if (n == 0) {
      lines <- c(vcf$left, read_lines(vcf$con, 0))
    } else {
      if (length(vcf$left) > n) {
        lines <- vcf$left[1:n]
        vcf$left <- vcf$left[(n+1):length(vcf$left)]
      } else if (length(vcf$left) == n) {
        lines <- vcf$left
        vcf$left <- vector(mode="character")
      } else {
        lines <- c(vcf$left, read_lines(vcf$con, n-length(vcf$left)))
        vcf$left <- vector(mode="character")
      }
    }
  } else {
    lines <- read_lines(vcf$con, n)
  }
  
  if(length(lines) == 0){
    return(data.table())
  }
  vcf$nvar <- vcf$nvar + length(lines)

  lines <- tstrsplit2(lines, split="\t", fixed=TRUE)
  if(length(lines) == 8 & length(vcf$header) == 9){
    vcf$header <- vcf$header[1:8]
  }
  if(length(vcf$header) != length(lines)){
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
  xs <- strsplit2(x,":",fixed=T)
  # ensure that the length of splited sample col is same as the splited format col.
  if(!all(sapply(xs, length) == FORMAT_length)){
    xs <- resize_list_string(xs, FORMAT_length, fill="")
  }
  unlist(xs)
}

str_to_ikv2 <- function(x, sep=";"){
  xs <- strsplit2(x, split=";", fixed=T)
  ikv <- data.table(i=rep(seq_along(x), sapply(xs, length)))
  kv <- strsplit2(unlist(xs, recursive=F, use.names=F), split="=", fixed=T)
  kv <- resize_list_string(kv, rep(2, length(kv)))
  ikv <- cbind(ikv, do.call(rbind, kv))
  setnames(ikv, c("i", "k", "v"))
  ikv
}

split_info <- function(x, info_keys=NULL, prefix="INFO."){
  info_ikv <- setDT(str_to_ikv(x, sep=";"))
  info_ikv <- info_ikv[!(i %in% which(x %in% c(".","",NA)))] # remove missing info
  
  if(!is.null(info_keys)){
    info_ikv <- info_ikv[k %in% info_keys]
  }
  if(nrow(info_ikv) == 0 ) {
    return(data.table())
  }

  setkey(info_ikv, i, k)
  info_ikv <- unique(info_ikv, by=c("i", "k"))

  info_dt <- dcast(info_ikv, i ~ k, value.var = "v")
  setnames(info_dt, c(".varid", colnames(info_dt)[-1]))
  idt <- data.table(.varid=1:length(x))
  info_dt <- info_dt[idt, on=".varid"][, -1]
  info_keys2 <- setdiff(info_keys, colnames(info_dt))
  if(length(info_keys2)>0){
    info_dt2 <- setDT(as.list(rep(NA_character_, length(info_keys2))))
    setnames(info_dt2, info_keys2)
    info_dt <- cbind(info_dt, info_dt2)
  }
  setnames(info_dt, paste0(prefix, colnames(info_dt)))
  info_dt
}

split_format <- function(x, format_keys=NULL, prefix="FORMAT."){
  stopifnot(inherits(x, "data.table"))

  if(!"FORMAT" %in% colnames(x)[1]){
    stop("not found FORMAT cols")
  }

  if(ncol(x)<2){
    stop(">=2 cols are needed")
  }

  # split sample format
  splited_format <- strsplit2(x$FORMAT, ":", fixed=TRUE)
  FORMAT_length <- sapply(splited_format, length)
  if(any(FORMAT_length == 0)){
    stop("Empty string is not allowed in the FORMAT col.")
  }
  format_ikv <- data.table(
    varid = rep(1:nrow(x), FORMAT_length),
    k = unlist(splited_format))

  format_ikv <- cbind(
    format_ikv,
    setDT(lapply(x[,-1], function(y){
      xs <- strsplit2(y, ":", fixed=T)
      xs <- resize_list_string(xs, FORMAT_length, fill="")
      unlist(xs)
    })))

  if(!is.null(format_keys)){
    format_ikv <-  format_ikv[k %in% format_keys]
  }
  setnames(format_ikv, c("varid", "k", colnames(x)[-1]))
  unique(format_ikv, by=c("varid","k"))
  format_dt <- dcast(format_ikv, varid ~ k, value.var = colnames(format_ikv)[3:ncol(format_ikv)])
  if(ncol(format_ikv) > 3){
    newnames <- matrix(stringr::str_match(colnames(format_dt)[-1], "(.+)_([^_]+)")[,-1], ncol=2)
    newnames <- paste0(newnames[,2], ".", newnames[,1])
    setnames(format_dt, c("varid", newnames))
  }else{
    setnames(format_dt, c("varid", paste0(colnames(format_dt)[-1], ".", colnames(format_ikv)[3])))
  }
  setkey(format_dt, varid)
  idt <- data.table(varid=1:nrow(x))
  format_dt <- format_dt[idt, on="varid"][, -1]

  format_keys2 <- setdiff(do.call(paste, c(as.list(expand.grid(format_keys, colnames(x)[-1])), sep=".")), colnames(format_dt))
  if(length(format_keys2)>0){
    format_dt[, get("format_keys2"):=NA_character_]
  }
  setnames(format_dt, paste0(prefix, colnames(format_dt)))
  format_dt
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
  if(ncol(lines) >= 9) {
    splited_format <- strsplit2(lines$FORMAT, ":", fixed=TRUE)
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

split_vars <- function(variants, info_keys=NULL, format_keys=NULL){
  dt <- data.table()
  if(length(variants)==0) return(dt)
  if(nrow(variants$C7) == 0 ) {
    return(dt)
  }else{
    dt <- variants$C7
  }

  if("INFO" %in% names(variants)){
    info_ikv <- variants$INFO
    if(!is.null(info_keys)){
      info_ikv <- info_ikv[k %in% info_keys]
    }
    info_dt <- dcast(info_ikv, varid ~ k, value.var = "v")
    setnames(info_dt, c("varid", paste0("INFO.", setdiff(colnames(info_dt),"varid"))))
    setkey(info_dt, varid)
    dt <- merge(dt, info_dt, by="varid", all.x=T, suffixes=c("", "."))
  }

  if("FORMAT" %in% names(variants)){
    format_ikv <- variants$FORMAT
    if(!is.null(format_keys)){
      format_ikv <-  format_ikv[k %in% format_keys]
    }
    unique(format_ikv, by=c("varid","k"))
    format_dt <- dcast(format_ikv, varid ~ k, value.var = setdiff(colnames(variants$FORMAT), c("varid","k")))
    newnames <- matrix(stringr::str_match(setdiff(colnames(format_dt), "varid"),"(.+)_([^_]+)")[,-1], ncol=2)
    newnames <- paste("FORMAT", newnames[,2], newnames[,1], sep=".")
    setnames(format_dt, c("varid", newnames))
    setkey(format_dt, varid)

    dt <- merge(dt, format_dt, by="varid", all.x=T, suffixes=c("", "."))
  }

  dt
}

#' Read n variants from a opened vcf file.
#' @param vcf A vcf object returned by open_vcf
#' @param n The number of lines to be read
#' @param info_keys a subset of INFO keys, missing for not spliting, NULL for all
#' @param format_keys a subset of FORMAT keys, missing for not spliting, NULL for all
#' @param sample_ids a subset of sample ids, missing for all
#' @details read_vars returns a data.table where the INFO and FORMAT are splited.
#' INFO feilds are named as INFO.key, INFO.AF for example.
#' FORMAT feilds are named as FORMAT.key.sampleid, FORMAT.GT.sample1 for example.
read_vars <- function(
  vcf, 
  n=0L, 
  info_keys, 
  format_keys, 
  sample_ids, 
  alt_break=FALSE, 
  info_keys_break=c("AC", "AF"), 
  format_keys_break=c("GT", "AD"))
{
  info_keys_break <- paste0("INFO.", info_keys_break)
  format_keys_break <- paste0("FORMAT.", format_keys_break)
  
  if(is.character(vcf)){
    vcf <- open_vcf(vcf)
  }
  stopifnot(inherits(vcf, "vcf"))
  dt <- data.table()
  lines_dt <- read_body(vcf, n=n)
  if("POS" %in% colnames(lines_dt)){
    lines_dt[, POS:=as.integer(POS)]
  }
  if(nrow(lines_dt) == 0) return(dt)
  
  if(!missing(sample_ids)){
    if(!is.null(sample_ids)){
      sample_ids <- na.omit(unique(sample_ids))
      id <- sample_ids %in% colnames(lines_dt)[-(1:9)]
      if(any(!id)){
        stop(paste0("unrecognized sample_ids: ",
                    paste(sample_ids[!id], collapse=","),
                    ".\n must be in the: ",
                    paste(colnames(lines_dt)[-(1:9)], collapse=", ")))
      }
      lines_dt <- lines_dt[, c(colnames(lines_dt)[1:9], sample_ids), with=F]
    }else{
      lines_dt <- lines_dt[, 1:8, with=F]
    }
  }
  
  dt <- lines_dt[, 1:7]
  alts <- strsplit2(dt$ALT, ",", fixed=TRUE)
  alts_count <- sapply(alts, length)
  alt_break <- alt_break & any(alts_count > 1)
  if(alt_break){
    nms <- colnames(dt)
    dt <- setDT(lapply(nms, function(x){
      if(x == "ALT"){
        res <- unlist(alts, use.names=F)
      }else{
        res <- rep(dt[[x]], alts_count)
      }
      res
    }))
    setnames(dt, nms)
  }
  
  if(!missing(info_keys)){
    info_dt <- split_info(lines_dt$INFO, info_keys=info_keys)
    if(nrow(info_dt)>0) {
      if(alt_break){
        nms <- colnames(info_dt)
        info_dt <- setDT(lapply(nms, function(x){
          if(x %in% info_keys_break) {
            res <- strsplit2(info_dt[[x]], ",", fixed=TRUE)
            if(any(sapply(res, length) != alts_count)){
              res <- unlist(resize_list_string(res, alts_count, fill=""), use.names=F)
            }else{
              res <- unlist(res, use.names=F)
            }
          }else{
            res <- rep(info_dt[[x]], alts_count)
          }
          res
        }))
        setnames(info_dt, nms)
      }
      dt <- cbind(dt, info_dt)
    }else{
      dt[, INFO:=""]
    }
  }else{
    if(alt_break){
      dt[, INFO:=rep(lines_dt$INFO, alts_count)]
    }else{
      dt <- cbind(dt, INFO=lines_dt$INFO)
    }
  }
  
  if(ncol(lines_dt)>9){
    fmt_dt <- lines_dt[, 9:ncol(lines_dt), with=FALSE]
    
    if(!missing(format_keys)){
      fmt_dt <- split_format(fmt_dt, format_keys=format_keys)
    }
    if(alt_break){
      nms <- colnames(fmt_dt)
      if(missing(format_keys)){
        fmt_dt <- setDT(lapply(nms, function(x){
          rep(fmt_dt[[x]], alts_count)
        }))
      }else{
        fmt_dt <- setDT(lapply(nms, function(x){
          xk <- regmatches(x, regexec("FORMAT\\.[^.]+",  x))[[1]]
          if(xk %in% format_keys_break){
            res <- strsplit2(fmt_dt[[x]], ",", fixed=TRUE)
            if(any(sapply(res, length) != alts_count)){
              res <- unlist(resize_list_string(res, alts_count, fill=""), use.names=F)
            }else{
              res <- unlist(res, use.names=F)
            }
          }else{
            res <- rep(fmt_dt[[x]], alts_count)
          }
          res
        }))
      }
      setnames(fmt_dt, nms)
    }
    dt <- cbind(dt, fmt_dt)
  }
  dt
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
    INFO <- uniq_char(do.call(paste, c(lapply(info_cols, function(x){
      k <- sub("^INFO\\.", "", x)
      v <- paste(k, DT[[x]], sep="=")
      v[is.na(DT[[x]])] <- ""
      v
    }),sep=";")), ";")
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

