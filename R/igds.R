table_region_lines <- function(lines, type_detect){
  if(type_detect){
    dt <- fread(paste0(paste(lines, collapse="\n"),'\n'), header=F, sep="\t")
  }else{
    dt <- setDT(tstrsplit2(lines, split="\t", fixed=TRUE))
  }
  dt[, V1:=as.character(V1)]
  suppressWarnings(dt[, V2:=as.integer(V2)])
  suppressWarnings(dt[, V3:=as.integer(V3)])
  dt
}

table_position_lines <- function(lines, type_detect){
  if(type_detect){
    dt <- fread(paste0(paste(lines, collapse="\n"),'\n'), header=F, sep="\t")
  }else{
    dt <- setDT(tstrsplit2(lines, split="\t", fixed=TRUE))
  }
  dt[, V1:=as.character(V1)]
  suppressWarnings(dt[, V2:=as.integer(V2)])
  dt
}

seq_group <- function(n, each){
  rep(1L:ceiling(n/each), each=each)[1L:n]
}

#' Get the index of chromosome intervals and reset the V2 (start position) and V3 (end position) of x as offset coordinates.
#' @param x a data.table of and other information.
#' @param binsize the size of each bin
#' @return a list including the index and x where the positions were modified.
bin_region <- function(x, binsize){
  x[, i:=.I]
  x[, bini:=seq_group(.N, binsize)]
  # star1 and end1 are used to recover the real coordinates using offset
  index <- x[, .(start=start[1],
                 end=max(end),
                 start_base=start[1],
                 end_base=end[1],
                 i=i[1],
                 isize=.N),
             by=c("chr","bini")]
  index[, bini:=NULL]
  # convert to offset, which will save large memory when ...
  x[, c("start", "end"):=
      .(c(0L, start[-1] - start[-.N]),
        c(0L, end[-1] - end[-.N])),
    by=c("chr","bini")]

  x[, c("i", "bini"):=NULL]
  list(index=index, data=x)
}

#' Get the index of variant postions and reset the V2 (variant position) of x as offset coordinates.
#' @param x a data.table of variants and other information.
#' @param binsize the size of each bin
#' @return a list including the index and x where the positions were modified.
bin_variant <- function(x, binsize){
  stopifnot(!missing(x))
  stopifnot(!missing(binsize))
  x[, i:=.I]
  x[, bini:=seq_group(.N, binsize)]
  index <- x[, .(start=pos[1],
                 end=pos[.N],
                 pos_base=pos[1],
                 i=i[1],
                 isize=.N),
             by=c("chr","bini")]
  index[, bini:=NULL]

  # convert to offset
  x[, c("pos") := .(c(0L, pos[-1] - pos[-.N])), by=c("chr", "bini")]

  x[, c("i", "bini"):=NULL]
  list(index=index, data=x)
}

#'
write_bin <- function(binlist, gds, dir="data", compress_args){
  datanode <- gdsfmt::index.gdsn(gds, dir, silent=T)
  if(is.null(datanode)) datanode <- gdsfmt::add.gdsn(gds, dir, val=list())

  for(chr in names(binlist)){
    bin <- binlist[[chr]]
    chrnode <- gdsfmt::index.gdsn(datanode, chr, silent=T)

    if(is.null(chrnode)) chrnode <- gdsfmt::add.gdsn(datanode, chr, val=data.table())

    for(col in colnames(bin$data)){
      colnode <- gdsfmt::index.gdsn(chrnode, col, silent=T)
      if(is.null(colnode)){
        if(col %in% names(compress_args)){
          args <- compress_args[[col]]
          if(!"compress" %in% names(compress_args)){
            args[["compress"]] <- 'LZ4_RA'
          }
          args[["node"]] <- chrnode
          args[["name"]] <- col
          args[["val"]] <- bin$data[[col]]
          args[["check"]] <- FALSE
          do.call(gdsfmt::add.gdsn, args)
        }else{
          gdsfmt::add.gdsn(chrnode, col, val=bin$data[[col]], compress='LZ4_RA',check=FALSE)
        }
      }else{
        gdsfmt::append.gdsn(colnode, val=bin$data[[col]],check=F)
      }
    }
  }
}

# i is the location of vector
update_index <- function(binlist, index_l,chrstatus){
  for(chr in names(binlist)){
    bin <- binlist[[chr]]
    if(chr %in% names(chrstatus)) bin$index[, i:=i+chrstatus[[..("chr")]]]

    index_l[[length(index_l)+1L]] <- bin$index
  }
  index_l
}

update_chrstatus <- function(x, chrstatus){
  for(chr in names(x)){
    if(! chr %in% names(chrstatus)){
      chrstatus[[chr]] <- as.integer(nrow(x[[chr]]))
    }else{
      chrstatus[[chr]] <- chrstatus[[chr]] + nrow(x[[chr]])
    }
  }
  chrstatus
}

parse_header <- function(reader){
  mystop <- function(msg){
    close_file(reader)
    stop(msg)
  }
  
  repeat{
    l1 <- read_lines(reader, 1L)
    if(!startsWith(l1, "##")) break
  }

  if(length(l1) == 0) {
    mystop("empty file!")
  }
  
  l1 <- sub("^#", "", l1)
  
  if(startsWith(l1, "chr\tstart\tend")){
    dbtype <- "region"
  }else if(startsWith(l1, "chr\tpos\tref\talt")){
    dbtype <- "variant"
  }else if(startsWith(l1, "chr\tpos")){
    dbtype <- "position"
  }else{
    mystop("unkown header")
  }
  
  list(header=strsplit2(l1, "\t")[[1]],
       dbtype=dbtype)
}

#' Indexing a genomic coordinate (1-based) file
#' @param file three types of files are allowed.
#' region file:
#' #chr start end gene score
#' chr1 1 999 A 1.1
#' chr2 9 222 B 2.2
#' variant file:
#' #chr pos ref alt gene score (1-based)
#' 1 22 A T gene1 1.1
#' 2 44 AT A gene2 22
#' position file:
#' #chr pos gene score (1-based)
#' 1 22 gene1 1.1
#' 2 44 gene2 22
#' 
#' pv=6MB/s

#' @param file Input file
#' @param dbname The gds file that will be created.
#' @param binsize How many lines stored in a bin
#' @param nbinread The number of bins to be read in one time
#' @param compress_args List, a list of compress arguments for compression.
#'
indexdb <- function(
  file, 
  dbname=paste0(file, ".igds"),
  compress_args=list(), 
  binsize=10000L,
  nbinread=1L,
  progress=FALSE){
  
  dbname.tmp <- paste0(dbname, ".indexdb.tmp")
  
  nlines <- binsize * nbinread
  stopifnot(file.exists(file))
  
  reader <- open_file(file)
  header <- parse_header(reader=reader)
  dbtype <- header$dbtype
  header <- header$header
  
  stopifnot(dbtype %in% c("region", "variant", "position"))
  index_cols <- dbtype_to_index_cols(dbtype)

  stopifnot(is.list(compress_args))

  if(any(names(compress_args) %in% index_cols)){
    close_file(reader)
    stop("cols in compress_args should not be the indexed cols")
  }
  
  if(any(!names(compress_args) %in% header)){
    close_file(reader)
    stop("cols in compress_args should be the header")
  }

  for(i in setdiff(header, index_cols)){
    if(is.null(compress_args[[i]])){
      compress_args[[i]] <- list(compress="LZ4_RA")
    }else{
      if(is.null(compress_args[[i]][["compress"]]))
        compress_args[[i]][["compress"]] <- "LZ4_RA"
    }
  }

  lines <- read_lines(reader=reader, n=nlines)
  if(length(lines)<1) {
    close_file(reader)
    stop("cannot find data lines!")
  }

  # make sure dbname is closed
  close_gdsfile(dbname.tmp)
  gds <- gdsfmt::createfn.gds(dbname.tmp)
  gdsfmt::add.gdsn(gds, "dbtype", val=dbtype)
  gdsfmt::add.gdsn(gds, "header", val=header)
  gdsfmt::add.gdsn(gds, "index_cols", val=index_cols)
  gdsfmt::add.gdsn(gds, "binsize", val=binsize)
  gdsfmt::add.gdsn(gds, "compress_args", val=compress_args)

  chrstatus <- vector("integer") # how many bases were stored for each chromosome.
  index_l <- list()
  total_records <- 0
  type_dectect <- TRUE # only detect the data types in the first round
  repeat{
    invisible(gc())
    if(length(lines) > 0){
      total_records <- total_records + length(lines)
      if(progress){
        message(total_records, " records were stored\r", appendLF=FALSE)
        flush.console()
      }
      lines_table <- switch(
        dbtype,
        region=table_region_lines(lines, type_dectect),
        position=table_position_lines(lines, type_dectect),
        variant=table_position_lines(lines, type_dectect))
      type_dectect <- FALSE
      
      #lines_table <- lines_table[, header_num, with=F]
      setnames(lines_table, header)
      setkeyv(lines_table, index_cols) # sorted is required for bin_region and bin_variant
      
      # split by chr
      uniq_chrs <- unique(lines_table$chr)
      if(length(uniq_chrs)==1){
        lines_table <- list(lines_table)
        names(lines_table) <- uniq_chrs
      }else{
        lines_table <- split(lines_table, by="chr")# split by chrom
      }
      
      binlist <- switch(
        dbtype,
        region=lapply(lines_table, bin_region, binsize=binsize),
        position=lapply(lines_table, bin_variant, binsize=binsize),
        variant=lapply(lines_table, bin_variant, binsize=binsize))
      
      index_l <- update_index(binlist, index_l=index_l, chrstatus=chrstatus)
      chrstatus <- update_chrstatus(lines_table, chrstatus=chrstatus)
      write_bin(binlist, gds=gds, dir="data", compress_args=compress_args)
    }
    
    if(length(lines <- read_lines(reader=reader, n=nlines))==0) {break}
  }
  
  gdsfmt::add.gdsn(gds, "total_records", val=total_records)

  index <- rbindlist(index_l)
  setkey(index, chr, start, end)
  indexnode <- gdsfmt::add.gdsn(gds, "index", val=index, compress="LZ4",check=F)
  gdsfmt::put.attr.gdsn(indexnode, name="sorted", val=c("chr","start","end"))

  close_file(reader)
  gdsfmt::closefn.gds(gds)
  gdsfmt::cleanup.gds(normalizePath(dbname.tmp), verbose=FALSE)
  file.rename(dbname.tmp, dbname)
  invisible()
}

locate_bin <- function(DT, index, gds){
  if(!"pos2" %in% colnames(DT)){
    DT[, pos2:=pos]
  }
  # fetch bin locations and the offset
  DTdx <- foverlaps(DT, index, by.x=c("chr","pos","pos2"), by.y=c("chr","start","end"), nomatch=0L)

  DTdx
}

#' Open a genome database
opendb <- function(file){
  stopifnot(inherits(file, "character"))
  file <- normalizePath(file, mustWork=T)
  gds <- gdsfmt::openfn.gds(file, allow.duplicate=TRUE,allow.fork=TRUE)
  header <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "header"))
  binsize <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "binsize"))
  dbtype <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "dbtype"))
  compress_args <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "compress_args"))
  
  indexnode <- gdsfmt::index.gdsn(gds, "index")
  index_attr <- gdsfmt::get.attr.gdsn(indexnode)
  index <- gdsfmt::read.gdsn(indexnode)
  if(is.null(index_attr[['sorted']])){
    gdsfmt::closefn.gds(gds)
    stop("not find the key of the index or index is not sorted!")
  }else{
    attr(index, 'sorted') <-  index_attr[['sorted']]
  }
  
  if(!all.equal(key(index), c("chr","start","end"))){
    gdsfmt::closefn.gds(gds)
    stop("the database is not indexed with chr, start and end")
  }
  
  datanode <- gdsfmt::index.gdsn(gds, "data")
  allchrnms <- gdsfmt::ls.gdsn(datanode)
  if(length(allchrnms)==0) {
    gdsfmt::closefn.gds(gds)
    stop("Can not find chromosomes in the database.")
  }
  chrnodes <- lapply(allchrnms, function(x) {gdsfmt::index.gdsn(datanode, path=x)})
  names(chrnodes) <- allchrnms
  
  db <- list(
    file=file,
    gds = gds,
    dbtype = dbtype,
    header = header,
    compress_args = compress_args,
    binsize = binsize,
    index = index,
    chrnodes = chrnodes
  )
  class(db) <- "igds"
  db
}

close.igds <- function(db){
  stopifnot(inherits(db, "igds"))
  gdsfmt::closefn.gds(db$gds)
}

fetchdb <- function(
  db, 
  chr, 
  pos, 
  pos2, 
  ref, 
  alt, 
  select, 
  ops=NULL,
  max_bin=2L,
  mc.cores=1,
  nomatch=NA,
  ol.minoverlap=1L,
  ol.type="any",
  ol.mult="all",
  warn=F,
  allops = list(
    first=function(x){x[1]},
    last=function(x){x[length(x)]},
    mean=mean,
    max=max,
    min=min,
    paste=function(x)paste(x, collapse=","))){
  
  if(!inherits(db, "igds")){
    stop("db is not a igds object, please use opendb first")
  }
  
  if(!(db$dbtype %in% c("variant", "region", "position"))){
    stop("dbtype is not in variant, region, position ")
  }
  
  stopifnot(!missing(chr))
  chr <- as.character(chr)
  stopifnot(!missing(pos))
  stopifnot(length(chr)==length(pos))
  query_variables <- c("chr", "pos")
  if(!missing(pos2)){
    stopifnot(length(chr)==length(pos2))
    query_variables <- c(query_variables, "pos2")
  }
  if(!missing(ref)){
    stopifnot(length(chr)==length(ref))
    query_variables <- c(query_variables, "ref")
  }
  if(!missing(alt)){
    stopifnot(length(chr)==length(alt))
    query_variables <- c(query_variables, "alt")
  }

  if(!missing(alt) & missing(ref)){
    stop("ref is not defined when alt is using")
  }
  if(!missing(pos2) & !missing(ref)){
    stop("ref and pos2 should not be used at the same time")
  }

  if(is.null(max_bin)) max_bin <- ceiling(100000L/db$binsize)
  
  if(missing(select)){
    select <- setdiff(db$header, query_variables)
  }

  if(!all(select %in% db$header)){
    stop(paste("these selected cols:",
               paste(select[!select %in% db$header], collapse=","),
               "must be in the header of database:", paste(db$header, collapse=', ')))
  }

  if(!is.null(ops)){
    if(!all(select %in% names(ops))){
      stop("if ops is used, please set ops for every select cols")
    }
    ops <- ops[select]
    
    opsfunc <- list()
    for(i in names(ops)){
      if(is.function(ops[[i]])){
        opsfunc[[i]] <- ops[[i]]
      }else if(inherits(ops[[i]], "character")){
        if(ops[[i]] %in% names(allops)){
          opsfunc[[i]] <- allops[[ops[[i]]]]
        }else{
          opsfunc[[i]] <- get(ops[[i]])
        }
      }else{
        stop("ops must be a function or character")
      }
    }
  }

  isinchrs <- chr %in% names(db$chrnodes)
  if(warn) {
    if(any(!isinchrs)) {
      warning(paste("There are chromosomes not in the database:",
                    paste(unique(chr[!isinchrs]), collapse=", ")))
    }
  }

  if(sum(isinchrs)==0) {
    if(warn){
      warning('All queried chromosomes are not found in the database, maybe need check their names.')
    }
    return(data.table())
  }
  DT <- setDT(lapply(query_variables, function(x) {get(x)[isinchrs]}))
  setnames(DT, query_variables)
  setkeyv(DT, query_variables)

  # variant / region to region ####
  if(db$dbtype == "region"){
    get_DTdxdb <- function(DTdx){
      bins <- unique(DTdx[, .(chr, i, isize, start_base, end_base)])
      f <- function(bi){
        chr_node <- db$chrnodes[[bins$chr[bi]]]
        start_bi <- cumsum(gdsfmt::read.gdsn(gdsfmt::index.gdsn(chr_node, "start"),
                                      start=bins$i[bi],
                                      count=bins$isize[bi])) + bins$start_base[bi]
        end_bi <- cumsum(gdsfmt::read.gdsn(gdsfmt::index.gdsn(chr_node, "end"),
                                    start=bins$i[bi],
                                    count=bins$isize[bi])) + bins$end_base[bi]
        chr_bi <- rep(bins$chr[bi], length(start_bi))
        select_bi <- lapply(select, function(col){
          gdsfmt::read.gdsn(gdsfmt::index.gdsn(chr_node, col),
                            start=bins$i[bi],
                            count=bins$isize[bi])
        })
        setDT(c(list(chr_bi), list(start_bi), list(end_bi), select_bi))
      }
      DTdb <- lapply(1:nrow(bins), f) %>% rbindlist

      setnames(DTdb, paste0("db.", c("chr","start", "end", select)))
      setkeyv(DTdb, c("db.chr", "db.start", "db.end"))

      DTdxdb <- foverlaps(DTdx,
                          DTdb,
                        by.x=c("chr","pos","pos2"),
                        by.y=c("db.chr","db.start","db.end"), 
                        nomatch=0,
                        minoverlap=ol.minoverlap,
                        type=ol.type,
                        mult=ol.mult)
      DTdxdb <- DTdxdb[, c(query_variables, paste0("db.", c("start", "end", select))), with=F]
      DTdxdb
    }
  }
  # region to variant ####
  # query: chr pos pos2; dbtype: variant
  if(db$dbtype %in% c("position", "variant") && !missing(pos2)){
    get_DTdxdb <- function(DTdx){
      bins <- unique(DTdx[, .(chr, i, isize, pos_base)])
      f <- function(bi){
        chr_node <- db$chrnodes[[bins$chr[bi]]]
        pos_bi <- cumsum(gdsfmt::read.gdsn(gdsfmt::index.gdsn(chr_node, "pos"),
                                    start=bins$i[bi],
                                    count=bins$isize[bi])) + bins$pos_base[bi]
        chr_bi <- rep(bins$chr[bi], length(pos_bi))
        select_bi <- lapply(select, function(col){
          gdsfmt::read.gdsn(gdsfmt::index.gdsn(chr_node, col),
                            start=bins$i[bi],
                            count=bins$isize[bi])
        })
        setDT(c(list(chr_bi),list(pos_bi),list(pos_bi), select_bi))
      }
      DTdb <- lapply(1:nrow(bins), f) %>% rbindlist
      setnames(DTdb, paste0("db.", c(c("chr","pos","pos2"), select)))
      setkeyv(DTdb, paste0("db.", c("chr","pos","pos2")))

      DTdxdb <- foverlaps(DTdx,
                          DTdb,
                          by.x=c("chr","pos","pos2"),
                          by.y=paste0("db.", c("chr","pos","pos2")), 
                          nomatch=0,
                          minoverlap=ol.minoverlap,
                          type=ol.type,
                          mult=ol.mult) %>%
        .[, c(query_variables, paste0("db.", c("pos", select))), with=F]
      DTdxdb
    }
  }

  # variant to variant ####
  # query: chr pos;         dbtype: variant
  # query: chr pos ref;     dbtype: variant
  # query: chr pos ref alt; dbtype: variant
  if(db$dbtype %in% c("position", "variant") && missing(pos2)){
    get_DTdxdb <- function(DTdx){
      bins <- unique(DTdx[, .(chr, i, isize, pos_base)])
      f <- function(bi){
        chr_node <- db$chrnodes[[bins$chr[bi]]]
        pos_bi <- cumsum(gdsfmt::read.gdsn(gdsfmt::index.gdsn(chr_node, "pos"),
                                    start=bins$i[bi],
                                    count=bins$isize[bi])) + bins$pos_base[bi]
        chr_bi <- rep(bins$chr[bi], length(pos_bi))
        select_bi <- lapply(c(setdiff(query_variables, c("chr","pos")), select), function(col){
          gdsfmt::read.gdsn(gdsfmt::index.gdsn(chr_node, col),
                            start=bins$i[bi],
                            count=bins$isize[bi])
        })
        setDT(c(list(chr_bi),list(pos_bi),select_bi))
      }
      DTdb <- lapply(1:nrow(bins), f) %>% rbindlist
      setnames(DTdb, paste0("db.", c(query_variables, select)))
      setkeyv(DTdb, paste0("db.", query_variables))

      DTdxdb <- merge(DTdx, DTdb, by.x=query_variables, by.y=paste0("db.", query_variables)) %>%
        .[, c(query_variables, paste0("db.", select)), with=F]
      DTdxdb
    }
  }

  DTdx <- locate_bin(DT=DT, index=db$index)
  if(nrow(DTdx)==0L) {
    return(data.table())
  }

  # split
  setkey(DTdx, chr, i)
  DTdx[, splitgrp:=c(1L, rep(0, .N-1)), by=c("chr","i")] %>%
    .[, splitgrp:=ceiling(cumsum(splitgrp)/max_bin)]
  DTdx_l <- split(DTdx, by=c("splitgrp"))

  if(mc.cores < 2){
    DTdxdb <- lapply(DTdx_l, get_DTdxdb)
  }else{
    DTdxdb <- parallel::mclapply(DTdx_l, get_DTdxdb, mc.cores=mc.cores)
  }
  DTdxdb <- rbindlist(DTdxdb)
  if(nrow(DTdxdb) == 0 ) {
    return(data.table())
  }
  setkeyv(DTdxdb, query_variables)
  
  if(!is.null(ops)){
    DTdxdb <- DTdxdb[, lapply(paste0("db.", ..("select")), function(x){
      opsfunc[[sub("^db.","",x)]](get(x))
    }), by=query_variables]
    setnames(DTdxdb, c(query_variables, paste0("db.", select)))
  }

  if(is.na(nomatch) & !is.null(ops)){
    DT0 <- setDT(lapply(query_variables, function(x) get(x)))
    setnames(DT0, query_variables)
    for(col in paste0("db.", select)){
      DT0[DTdxdb, (get("col")):=get(..("col")), on=query_variables]
    }
    DTdxdb <- DT0
  }

  DTdxdb
}

#' close gds file handle or close all gds file hanldes connected to the file if a file name is provided.
#' @param x a gds object or a file
close_gdsfile <- function(x){
  if(inherits(x, "gds.class")){
    gdsfmt::closefn.gds(x)
  }
  if(inherits(x, "character")){
    rv <- .Call(gdsfmt:::gdsGetConnection)
    if (length(rv) > 0L) {
      for (i in seq_along(rv)) {
        names(rv[[i]]) <- c("filename", "id", "root", "readonly")
        class(rv[[i]]) <- "gds.class"
      }
    }
    else rv <- NULL
    
    for(i in rv){
      if(normalizePath(i$filename) == normalizePath(x)){
        gdsfmt::closefn.gds(i)
      }
    }
  }
  invisible()
}

dbtype_to_index_cols <- function(dbtype){
  switch(
    dbtype,
    position=c("chr","pos"),
    region=c("chr","start","end"),
    variant=c("chr","pos","ref","alt")
  )
}

print.igds <- function(db){
  stopifnot(inherits(db, "igds"))
  
  cat(paste0("file: ", db$file, "\n"))
  cat(paste0("dbtype: ", db$dbtype, "\n"))
  index_cols <- dbtype_to_index_cols(db$dbtype)
  cat(paste0("index cols: ", paste0(index_cols, collapse=", "), "\n"))
  cat(paste0("feature cols: ", paste(setdiff(db$header, index_cols), collapse=", "), "\n"))
  cat(paste0("total bins: ", nrow(db$index), "\n"))
  cat(paste0("bin size: ", nrow(db$binsize), "\n"))
  cat(paste0("total records: ", sum(as.numeric(db$index$isize)), "\n"))
}

head.igds <- function(db, n=10){
  stopifnot(inherits(db, "igds"))
  
  a <- fetchdb(db, chr=db$index[["chr"]][[1]], 
               pos=db$index[["start"]][[1]],
               pos2=db$index[["end"]][[1]])
  a <- a[, -c("pos", "pos2")]
  
  setnames(a, sub("db.", "", colnames(a)))
  head(a, n=n)
}
