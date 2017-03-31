get_filecon <- function(file){
  if(grepl(".gz$", file)){
    con <- gzfile(file, open="rb")
  }else{
    con <- file(file, open="rb")
  }
  con
}

open_vcf <- function(file, header_max=1000){
  #' read vcf file as a list cluding header
  #' @param file path of vcf file
  #' @param header_max the max lines to search for the vcf header, if 0, skips reading header, used only for no header lines
  #'

  vcf <- new.env(parent = emptyenv())
  vcf$filepath <- normalizePath(file, mustWork=T)
  
  # read header
  con <- get_filecon(file)
  suppressWarnings(header <- readLines(con, n=header_max))
  close(con)
  
  header_pos <- which(substr(header, 1,6)=="#CHROM")
  if(length(header_pos)==0){
    stop(paste("can not find header in the first", header_max,  "lines"))
  }else{
    vcf$header <- header[1:(header_pos)]
    vcf$header_names <- strsplit(sub("^#","",header[header_pos]), "\t")[[1]]
  }

  if(length(vcf$header_names) < 8){
    stop("vcf cols is not < 8")
  }
  if(!all.equal(vcf$header_names[1:8] , c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"))){
    stop("incorrect header names")
  }
  
  # get file connection skipping header
  vcf$filecon <- get_filecon(file)
  suppressWarnings(a <- readLines(vcf$filecon, n=length(vcf$header))) # skip header
  
  # init variant count
  vcf$variant_count <- 0L
  
  class(vcf) <- c("vcf", "environment")
  vcf
}

close_vcf <- function(vcf){
  close(vcf$filecon)
}

read_vcf <- function(file, get_info=TRUE){
  #' read vcf file variant lines (data.table)
  #'
  #' @param vcf an env obj
  
  vcf <- open_vcf(file)
  close(vcf$filecon)
  
  which_cat <- "cat"
  if(grepl(".gz$", vcf$filepath))
    which_cat <-  "zcat"
  
  lines <- fread(paste(which_cat, vcf$filepath), sep="\t", colClasses=c("char"),
                 skip=length(vcf$header), header=F)
  if(ncol(lines) < 8)
    stop("col num is < 8")
  vcf$header_names <- vcf$header_names[1:ncol(lines)]
  colnames(lines) <- vcf$header_names
  
  if(nrow(lines) == 0)
    stop("can not find variants")
  
  lines[, POS:=as.integer(POS)]
  
  vcf$lines <- lines
  
  if(get_info)
    parse_info(vcf)
  vcf
}

read_chunk <- function(vcf, chunk_size=1e5){
  #' read vcf file variant lines (data.table)
  #'
  #' @param vcf an env obj
  #' 

  vcf$lines <- NULL # empty old lines
  
  suppressWarnings(lines <- readLines(vcf$filecon, n=chunk_size))
  if(length(lines) == 0){
    return(FALSE)
  }
  lines <- strsplit_wide(lines, "\t", fixed=TRUE)
  
  if(ncol(lines) < 8)
    stop("col num is < 8")
  
  stopifnot(ncol(lines) <= length(vcf$header_names)) # supporting reductant header names
  
  setnames(lines,  vcf$header_names[1:ncol(lines)])
  
  if(nrow(lines) == 0)
  {
    print("there are not variants anymore!")
    close(vcf$filecon)
    return(NULL)
  }
  
  lines[, POS:=as.integer(POS)]
  
  vcf$lines <- lines
  
  parse_info(vcf)
  
  vcf$variant_count = vcf$variant_count + nrow(lines)
  
  return(TRUE)
}

write_vcf <- function(vcf, file, append=FALSE){
  if(!append & !is.null(v?svcf$header))
    write(vcf$header, file)
  fwrite(vcf$lines, file, append=T, sep="\t")
}

parse_info <- function(vcf){
  #' convert key-value sets info into a data.table
  a <- strsplit_long(vcf$lines$INFO, ";", fixed=TRUE)
  setnames(a, c("varid","kv"))
  a[, c("k", "v"):=tstrsplit(kv, "=", type.convert=F, fixed=T, fill="")]
  a[, kv:=NULL]
  setkey(a, k, varid)
  vcf$INFO <- a
}

parse_sample_format <- function(vcf){
  #' extract sample format
  #'
  #' @param vcf vcf env obj

  if(ncol(vcf$lines) >= 10){
    if("FORMAT" != colnames(vcf$lines)[9]){
      stop("col 9 is not FORMAT")
    }
    # FORMAT <- strsplit_long(vcf$lines$FORMAT, ":", fixed=T)
    # setnames(FORMAT, c("varid", "k"))
    FORMAT <- data.table(
      varid=rep(1:nrow(vcf$lines), each=4),
      k=unlist(transpose(tstrsplit(vcf$lines$FORMAT, ":", fixed=TRUE, keep=1:4))))
    sample_names <- colnames(vcf$lines)[-(1:9)]
    #FORMAT[, (sample_names):=lapply(sample_names, function(x) unlist(strsplit(vcf$lines[[x]], ":", fixed=T)))]
    FORMAT[, (sample_names):=lapply(sample_names, function(x) unlist(transpose(tstrsplit(vcf$lines[[x]], ":", fixed=TRUE, keep=1:4))))]
    setkey(FORMAT, k, varid)

    vcf$FORMAT <- FORMAT
  }
}

get_genotype <- function(vcf){
  ### get genotypes
  
  if(is.null(vcf$FORMAT))
    parse_sample_format(vcf)
  
  gt <-  vcf$FORMAT["GT", on="k"]
  if(nrow(gt) > 0){
    gt <- gt[, -"k", with=F]
    sample_names <- setdiff(colnames(gt), "varid")
    gt[, (sample_names):=lapply(.SD[,sample_names,with=F], function(x) {
      x <- sub("|", "/", x, fixed=TRUE)
      a <- tstrsplit(x,"/", fixed=TRUE)
      y <- rep(0L, length(x))
      y[a[[1]] != a[[2]]] <- 1L
      y[a[[1]] == a[[2]] & a[[1]] >0] <- 2L
      y[a[[1]] == "." | is.na(a[[2]]) | a[[2]] == "."] <- -9L # missing values are set to -9L
      y
    })]
    vcf$genotype <- gt
  }
  
  invisible(gt)
}


