table_region_lines <- function(lines){
  dt <- fread(paste0(paste(lines, collapse="\n"),'\n'), header=F, sep="\t")
  dt[, V1:=as.character(V1)]
  dt[, V2:=as.integer(V2)]
  dt[, V3:=as.integer(V3)]
  dt
}

table_position_lines <- function(lines){
  dt <- fread(paste0(paste(lines, collapse="\n"),'\n'), header=F, sep="\t")
  dt[, V1:=as.character(V1)]
  dt[, V2:=as.integer(V2)]
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

parse_header <- function(reader, index_cols){
  n <- length(index_cols)
  repeat{
    l1 <- read_lines(reader, 1L)
    if(!startsWith(l1, "#")) break
  }

  if(length(l1) == 0) {
    close_file(reader)
    stop("empty file!")
  }

  if(!startsWith(l1, paste(index_cols, collapse="\t"))){
    close_file(reader)
    stop(paste("the leading cols must be:",
               paste(index_cols, collapse="\t")))
  }

  l1_split <- strsplit(l1, "\t")[[1]]

  if(length(l1_split) <= n) {
    close_file(reader)
    stop(paste("the first line of the file must be tab-separated and >=",
               n + 1, "cols are needed!"))
  }

  l1_split
}

#' create indexed gds database of a region or variant file
#' @param file a region or variant file (1-based)
#' region file like this:
#' #chr start end gene score
#' chr1 1 999 A 1.1
#' chr2 9 222 B 2.2
#' variant file should like this:
#' #chr pos ref alt gene score (1-based)
#' 1 22 A T gene1 1.1
#' 2 44 AT A gene2 22
#' @param dbtype one of c("region", "variant")
#' @param gdsfile the gds file that will be created.
#' @param binsize How many lines stored in a bin
#' @param nbinread
#' @param compress_args
#' @param storage
#'
indexdb <- function(file, dbtype, gdsfile=paste0(file, ".gds"),
                    compress_args, binsize=10000L, nbinread=1L){
  nlines <- binsize * nbinread
  stopifnot(file.exists(file))
  stopifnot(!missing(dbtype))
  stopifnot()
  stopifnot(dbtype %in% c("region", "variant", "position"))
  index_cols <- switch(
    dbtype,
    region=c("chr","start","end"),
    position=c("chr","pos"),
    variant=c("chr","pos","ref","alt"))

  reader <- open_file(file)
  header <- parse_header(reader=reader, index_cols=index_cols)

  if(missing(compress_args)){
    compress_args <- list()
    for(i in setdiff(header, index_cols)){
      compress_args[[i]] <- list(compress="LZ4_RA")
    }
  }
  
  stopifnot(is.list(compress_args))

  if(any(names(compress_args) %in% index_cols)){
    close_file(reader)
    stop("cols in compress_args should not be the indexed cols")
  }
  
  if(any(!names(compress_args) %in% header)){
    close_file(reader)
    stop("cols in compress_args should be the header")
  }
  header_num <- match(c(index_cols, names(compress_args)), header)
  header <- c(index_cols, names(compress_args))

  for(i in names(compress_args)){
    if(is.null(compress_args[[i]][["compress"]]))
      compress_args[[i]][["compress"]] <- "LZ4_RA"
  }

  lines <- read_lines(reader=reader, n=nlines)
  if(length(lines)<1) {
    close_file(reader)
    stop("cannot find data lines!")
  }

  # make sure gdsfile is closed
  close_gdsfile(gdsfile)
  gds <- gdsfmt::createfn.gds(gdsfile)
  gdsfmt::add.gdsn(gds, "dbtype", val=dbtype)
  gdsfmt::add.gdsn(gds, "header", val=header)
  gdsfmt::add.gdsn(gds, "index_cols", val=index_cols)
  gdsfmt::add.gdsn(gds, "binsize", val=binsize)
  gdsfmt::add.gdsn(gds, "compress_args", val=compress_args)

  chrstatus <- vector("integer") # how many bases were stored for each chromosome.
  index_l <- list()
  repeat{
    lines_table <- switch(
      dbtype,
      region=table_region_lines(lines),
      position=table_position_lines(lines),
      variant=table_position_lines(lines))

    lines_table <- lines_table[, header_num, with=F]
    setnames(lines_table, header)
    setkeyv(lines_table, index_cols) # sorted is required for bin_region and bin_variant

    lines_table <- split(lines_table, by="chr")# split by chrom

    binlist <- switch(
      dbtype,
      region=lapply(lines_table, bin_region, binsize=binsize),
      position=lapply(lines_table, bin_variant, binsize=binsize),
      variant=lapply(lines_table, bin_variant, binsize=binsize))

    index_l <- update_index(binlist, index_l=index_l, chrstatus=chrstatus)
    chrstatus <- update_chrstatus(lines_table, chrstatus=chrstatus)
    write_bin(binlist, gds=gds, dir="data", compress_args=compress_args)

    if(length(lines <- read_lines(reader=reader, n=nlines))==0) {break}
  }

  index <- rbindlist(index_l)
  setkey(index, chr, start, end)
  indexnode <- gdsfmt::add.gdsn(gds, "index", val=index, compress="LZ4",check=F)
  gdsfmt::put.attr.gdsn(indexnode, name="sorted", val=c("chr","start","end"))

  close_file(reader)
  gdsfmt::closefn.gds(gds)
  gdsfmt::cleanup.gds(normalizePath(gdsfile), verbose=FALSE)
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

fetchdb <- function(
  gds, chr, pos, pos2, ref, alt, select, ops=NULL,
  max_bin=2L, mc.cores=1,
  nomatch=NA,
  warn=F,
  allops = list(
    first=function(x){x[1]},
    last=function(x){x[length(x)]},
    mean=mean,
    max=max,
    min=min,
    paste=function(x)paste0(x, collapse=",")))
{
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

  closegds <- F
  if(inherits(gds, "character")){
    gds <- gdsfmt::openfn.gds(gds, allow.duplicate=TRUE,allow.fork=TRUE)
    closegds <- T
  }

  header <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "header"))
  if(missing(select)){
    select <- setdiff(header, query_variables)
  }

  if(!all(select %in% header)){
    if(closegds) gdsfmt::closefn.gds(gds)
    stop(paste("these selected cols:",
               paste(select[!select %in% header], collapse=","),
               "must be in the header of database:", paste(header, collapse=', ')))
  }

  if(!is.null(ops)){
    if(length(select)!=length(ops)){
      if(closegds) gdsfmt::closefn.gds(gds)
      stop("cols must be equal with select")
    }

    if(!all(ops %in% names(allops))){
      if(closegds) gdsfmt::closefn.gds(gds)
      stop(paste("ops must be in:", paste(names(allops), collapse=",")))
    }

    opsfunc <- allops[ops]
    names(opsfunc) <- select
  }

  binsize <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "binsize"))
  if(is.null(max_bin)) max_bin <- ceiling(100000L/binsize)

  indexnode <- gdsfmt::index.gdsn(gds, "index")
  index_attr <- gdsfmt::get.attr.gdsn(indexnode)
  index <- gdsfmt::read.gdsn(indexnode)
  if(is.null(index_attr[['sorted']])){
    if(closegds) gdsfmt::closefn.gds(gds)
    stop("not find the key of the index!")
  }else{
    attr(index, 'sorted') <-  index_attr[['sorted']]
  }
  if(!all.equal(key(index), c("chr","start","end"))){
    if(closegds) gdsfmt::closefn.gds(gds)
    stop("the database is not indexed with chr, start and end")
  }

  datanode <- gdsfmt::index.gdsn(gds, "data")
  allchrnms <- gdsfmt::ls.gdsn(datanode)
  if(length(allchrnms)==0) {
    if(closegds) gdsfmt::closefn.gds(gds)
    stop("Can not find chromosomes in the database.")
  }
  chrnodes <- lapply(allchrnms, function(x) {gdsfmt::index.gdsn(datanode, path=x)})
  names(chrnodes) <- allchrnms

  isinchrs <- chr %in% allchrnms
  if(warn) {
    if(any(!isinchrs)) {
      warning(paste("There are chromosomes not in the database:",
                    paste(unique(chr[!isinchrs]), collapse=", ")))
    }
  }

  if(sum(isinchrs)==0) {
    if(closegds) gdsfmt::closefn.gds(gds)
    if(warn){
      warning('All queried chromosomes are not found in the database, maybe need check their names.')
    }
    return(data.table())
  }
  DT <- setDT(lapply(query_variables, function(x) {get(x)[isinchrs]}))
  setnames(DT, query_variables)
  setkeyv(DT, query_variables)

  dbtype <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "dbtype"))

  # variant / region to region ####
  if(dbtype == "region"){
    get_DTdxdb <- function(DTdx){
      bins <- unique(DTdx[, .(chr, i, isize, start_base, end_base)])
      f <- function(bi){
        chr_node <- chrnodes[[bins$chr[bi]]]
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
      db <- lapply(1:nrow(bins), f) %>% rbindlist

      setnames(db, paste0("db.", c("chr","start", "end", select)))
      setkeyv(db, c("db.chr", "db.start", "db.end"))

      DTdxdb <- foverlaps(DTdx,
                        db,
                        by.x=c("chr","pos","pos2"),
                        by.y=c("db.chr","db.start","db.end"), nomatch=0)
      DTdxdb <- DTdxdb[, c(query_variables, paste0("db.", c("start", "end", select))), with=F]
      DTdxdb
    }
  }
  # region to variant ####
  # query: chr pos pos2; dbtype: variant
  if(dbtype %in% c("position", "variant") && !missing(pos2)){
    get_DTdxdb <- function(DTdx){
      bins <- unique(DTdx[, .(chr, i, isize, pos_base)])
      f <- function(bi){
        chr_node <- chrnodes[[bins$chr[bi]]]
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
      db <- lapply(1:nrow(bins), f) %>% rbindlist
      setnames(db, paste0("db.", c(c("chr","pos","pos2"), select)))
      setkeyv(db, paste0("db.", c("chr","pos","pos2")))

      DTdxdb <- foverlaps(DTdx,
                          db,
                          by.x=c("chr","pos","pos2"),
                          by.y=paste0("db.", c("chr","pos","pos2")), nomatch=0) %>%
        .[, c(query_variables, paste0("db.", c("pos", select))), with=F]
      DTdxdb
    }
  }

  # variant to variant ####
  # query: chr pos;         dbtype: variant
  # query: chr pos ref;     dbtype: variant
  # query: chr pos ref alt; dbtype: variant
  if(dbtype == c("position", "variant") && missing(pos2)){
    get_DTdxdb <- function(DTdx){
      bins <- unique(DTdx[, .(chr, i, isize, pos_base)])
      f <- function(bi){
        chr_node <- chrnodes[[bins$chr[bi]]]
        pos_bi <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(chr_node, "pos"),
                                    start=bins$i[bi],
                                    count=bins$isize[bi]) + bins$pos_base[bi]
        chr_bi <- rep(bins$chr[bi], length(pos_bi))
        select_bi <- lapply(c(setdiff(query_variables, c("chr","pos")), select), function(col){
          gdsfmt::read.gdsn(gdsfmt::index.gdsn(chr_node, col),
                            start=bins$i[bi],
                            count=bins$isize[bi])
        })
        setDT(c(list(chr_bi),list(pos_bi),select_bi))
      }
      db <- lapply(1:nrow(bins), f) %>% rbindlist
      setnames(db, paste0("db.", c(query_variables, select)))
      setkeyv(db, paste0("db.", query_variables))

      DTdxdb <- merge(DTdx, db, by.x=query_variables, by.y=paste0("db.", query_variables)) %>%
        .[, c(query_variables, paste0("db.", select)), with=F]
      DTdxdb
    }
  }

  DTdx <- locate_bin(DT=DT, index=index)
  if(nrow(DTdx)==0L) {
    if(closegds) gdsfmt::closefn.gds(gds)
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
    if(closegds) gdsfmt::closefn.gds(gds)
    return(data.table())
  }
  setkeyv(DTdxdb, query_variables)

  if(!is.null(ops)){
    DTops <- DTdxdb[, lapply(paste0("db.", select), function(x){
      opsfunc[[x]](get(x))
    }), by=query_variables]
    setnames(DTops, c(query_variables, paste0("db.", select)))

    if(is.na(nomatch)){
      DT0 <- setDT(lapply(query_variables, get))
      setnames(DT0, query_variables)
      for(col in paste0("db.", select)){
        DT0[DTops, (col):=get(col), on=query_variables]
      }
      DTdxdb <- DT0
    }
  }

  if(closegds) gdsfmt::closefn.gds(gds)
  DTdxdb
}

function (closeall = FALSE, verbose = TRUE)
{
  stopifnot(is.logical(closeall))
  stopifnot(is.logical(verbose))
  rv <- .Call(gdsGetConnection)
  if (length(rv) > 0L) {
    nm <- NULL
    rd <- NULL
    for (i in seq_along(rv)) {
      names(rv[[i]]) <- c("filename", "id", "root", "readonly")
      class(rv[[i]]) <- "gds.class"
      nm <- c(nm, rv[[i]]$filename)
      rd <- c(rd, rv[[i]]$readonly)
    }
    if (verbose & !closeall) {
      print(data.frame(FileName = nm, ReadOnly = rd, State = rep("open",
                                                                 length(rd))))
    }
  }
  else rv <- NULL
  if (closeall & !is.null(rv)) {
    if (verbose) {
      print(data.frame(FileName = nm, ReadOnly = rd, State = rep("closed",
                                                                 length(rd))))
    }
    for (i in seq_along(rv)) closefn.gds(rv[[i]])
    rv <- NULL
  }
  invisible(rv)
}

#' close gds file by file path or a object of gds class,
#' if x is a file path,  all file handles connected to it will be closed.
#' @param x a gds object or a file path
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
