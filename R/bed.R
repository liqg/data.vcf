lines2bed <- function(lines){
    bed <- fread(paste(lines, collapse="\n"), header=F, sep="\t")
    bed[, V1:=as.character(V1)]
    bed[, V2:=as.integer(V2)+1L] # covert to 1-based
    bed[, V3:=as.integer(V3)]
    bed
}

seqgroup <- function(each, len){
    rep(1L:ceiling(len/each), each=each)[1L:len]
}

bedindex <- function(bed, binsize){
    bed[, i:=.I]
    bed[, bini:=seqgroup(binsize, .N)]
    index <- bed[,.(start=V2[1],end=V3[1], end2=V3[.N],i=i[1]),by=c("V1","bini")]
    index[, bini:=NULL]

    # convert to offset
    bed[, c("V2","V3"):=
            .(cumsum(c(0L, V2[-1] - V2[-.N])),
              cumsum(c(0L, V3[-1] - V3[-.N]))),
        by=c("V1","bini") ]
    
    bed[, c("i", "bini"):=NULL]
    list(index=index, bed=bed)
}

writetogds <- function(binlist, gds, dir="data", storage, compress){
    datanode <- gdsfmt::index.gdsn(gds, dir, silent=T)
    if(is.null(datanode)) datanode <- gdsfmt::add.gdsn(gds, dir, val=list())
    
    for(chr in names(binlist)){
        bin <- binlist[[chr]]
        chrnode <- gdsfmt::index.gdsn(datanode, chr, silent=T)

        if(is.null(chrnode)) chrnode <- gdsfmt::add.gdsn(datanode, chr, val=data.table())
        
        for(col in colnames(bin$bed)){
            colnode <- gdsfmt::index.gdsn(chrnode, col, silent=T)
            if(is.null(colnode)){
                if(!is.null(storage) & col %in% names(storage)){
                    gdsfmt::add.gdsn(chrnode, col, val=bin$bed[[col]], compress=compress, storage=storage[[col]], check=FALSE)
                }else{
                    gdsfmt::add.gdsn(chrnode, col, val=bin$bed[[col]], compress=compress,check=FALSE)
                }
            }else{
                gdsfmt::append.gdsn(colnode, val=bin$bed[[col]],check=F)
            }
        }
    }
}

update_index <- function(binlist, index_l,chrstatus){
    for(chr in names(binlist)){
        bin <- binlist[[chr]]
        if(chr %in% names(chrstatus)) bin$index[, i:=i+chrstatus[[chr]]]
        
        index_l[[length(index_l)+1L]] <- bin$index
    }
    index_l
}

update_chrstatus <- function(bedlist, chrstatus){
    for(chr in names(bedlist)){
        if(! chr %in% names(chrstatus)){
            chrstatus[[chr]] <- as.integer(nrow(bedlist[[chr]]))
        }else{
            chrstatus[[chr]] <- chrstatus[[chr]] + nrow(bedlist[[chr]])
        }
    }
    chrstatus
}

bed2gds <- function(bedfile, gdsfile, binsize=10000L, nbinread=1L, compress="LZ4_RA", storage=NULL){
    # stopifnot(file.exists(bedfile))
    nlines <- binsize*nbinread
    
    reader <- open_file(bedfile)
    l1 <- read_lines(reader,1)
    if(length(l1)!=1) stop("empty bed file!")
    l1_split <- strsplit(l1, "\t")[[1]]
    if(length(l1_split) < 4) stop(">= 4 cols are columns are needed!")
    if(substr(l1,1,1)=="#"){
        header <- l1_split
        header[1] <- sub("#", "", header[1])
        withheader <- T
    }else{
        header <- paste0("V", 1:length(l1_split))
        withheader <- F
    }
    
    if(withheader){
        lines <- read_lines(reader=reader, n=nlines)
    }else{
        lines <- c(l1,read_lines(reader=reader, n=nlines-1))
    }
    
    if(length(lines)<1) stop("cannot find data!")
    
    gds <- gdsfmt::createfn.gds(gdsfile)
    gdsfmt::add.gdsn(gds, "header", val=header)
    gdsfmt::add.gdsn(gds, "binsize", val=binsize)

    chrstatus <- vector("integer")
    index_l <- list()
    
    bedlist <- split(lines2bed(lines), by="V1")# split by chrom
    binlist <- lapply(bedlist, bedindex, binsize=binsize)
    index_l <- update_index(binlist, index_l=index_l, chrstatus=chrstatus)
    chrstatus <- update_chrstatus(bedlist, chrstatus=chrstatus)
    writetogds(binlist, gds=gds, dir="data", storage=storage, compress=compress)

    while(length(lines <- read_lines(reader=reader, n=nlines))>0){
        bedlist <- split(lines2bed(lines), by="V1")# split by chrom
        binlist <- lapply(bedlist, bedindex, binsize=binsize)
        index_l <- update_index(binlist, index_l=index_l, chrstatus=chrstatus)
        chrstatus <- update_chrstatus(bedlist, chrstatus=chrstatus)
        writetogds(binlist, gds=gds, dir="data", storage=storage, compress=compress)
    }
    
    index <- rbindlist(index_l)
    setkey(index, V1, start, end2)
    indexnode <- gdsfmt::add.gdsn(gds, "index", val=index, compress=compress,check=F)
    gdsfmt::put.attr.gdsn(indexnode, name="sorted", val=c("V1","start","end2"))
    
    close_file(reader)
    gdsfmt::closefn.gds(gds)
    gdsfmt::cleanup.gds(normalizePath(gdsfile))
    invisible()
}

bedanno <- function(chr, pos, gds, cols, newnames=paste0("V", cols), ops=NULL, 
                    max.bin=NULL, mc.cores=1,
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
    if(length(cols)!=length(newnames)){
        stop("cols must be equal with ops")
    }
    Vcols <- paste0("V", cols)
    
    if(is.null(ops)) {
        ops <- rep("first", length(cols))
    }
    if(length(cols)!=length(ops)){
        stop("cols must be equal with ops")
    }

    if(!all(ops %in% names(allops))){
        stop(paste("ops must be in:", paste(names(allops), collapse=",")))
    }
    
    opsfunc <- allops[ops]
    names(opsfunc) <- Vcols
    
    closegds <- F
    if(inherits(gds, "character")){
        gds <- gdsfmt::openfn.gds(gds, allow.duplicate=TRUE,allow.fork=TRUE)
        closegds <- T
    }   
    datanode <- gdsfmt::index.gdsn(gds, "data")
    binsize <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "binsize"))
    if(is.null(max.bin)) max.bin <- ceiling(100000L/binsize)
    index <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "index"))
    if(is.null(key(index))) setkey(index, V1, start, end2)
    if(!all.equal(key(index), c("V1","start","end2"))){
        setkey(index, V1, start, end2)
    }
    allchrnms <- gdsfmt::ls.gdsn(datanode)
    if(length(allchrnms)==0) stop("Can not find chromosomes in the database.")
    chrnodes <- lapply(allchrnms, function(x) {gdsfmt::index.gdsn(datanode, path=x)})
    names(chrnodes) <- allchrnms
    isinchrs <- chr %in% allchrnms
    if(warn) {
        if(any(!isinchrs)) warning("There are chromosomes being not in the database!")
    }
    
    DT <- data.table(chr=chr[isinchrs], pos=pos[isinchrs])
    if(nrow(DT) == 0) return(DT)
    
    setkey(DT, chr, pos)
    DT <- unique(DT)
    
    DT[, pos2:=pos]
    DTdx <- foverlaps(DT, index, by.x=c("chr","pos","pos2"), by.y=c("V1","start","end2"), nomatch=0L) %>%
        .[, V2:=pos - start] %>%
        .[, V3:=pos - end] %>%
        .[,.(chr, pos, i, V2, V3)]
    
    if(nrow(DTdx)==0L) return(data.table())
    
    get_DTdxdb <- function(DTdx){
        bins <- unique(DTdx[, .(chr, i)])
        db <- lapply(1:nrow(bins), function(bi){
            chrnode <- chrnodes[[bins$chr[bi]]]
            setDT(lapply(paste0("V", c(1:3, cols)), function(col){
                gdsfmt::read.gdsn(gdsfmt::index.gdsn(chrnode, col), start=1L, count=binsize)
            }))}) %>% rbindlist
        
        setnames(db, paste0("V", c(1:3, cols)))
        setkey(db, V1,V2,V3)
            
        foverlaps(DTdx[,.(chr, V2,V3,pos)],
                  db,
                  by.x=c("chr","V2","V3"),
                  by.y=c("V1","V2","V3"), nomatch=0)[, c("chr","pos", Vcols), with=F]
    }
    
    setkey(DTdx, chr, i)
    DTdx[, splitgrp:=c(1L, rep(0, .N-1)), by=c("chr","i")] %>%
        .[, splitgrp:=ceiling(cumsum(splitgrp)/max.bin)]
    DTdx_l <- split(DTdx, by=c("splitgrp"))
    if(mc.cores < 2){
        DTdxdb <- lapply(DTdx_l, get_DTdxdb)
    }else{
        DTdxdb <- parallel::mclapply(DTdx_l, get_DTdxdb, mc.cores=mc.cores)
    }
    DTdxdb <- rbindlist(DTdxdb)
    
    if(nrow(DTdxdb)>0){
        setkey(DTdxdb, chr, pos)
        DTops <- DTdxdb[, lapply(Vcols, function(x){
            opsfunc[[x]](get(x))
        }) ,by=c("chr","pos")]
        setnames(DTops, c("chr","pos", newnames))
        
        if(!is.na(nomatch) && nomatch==0){
          return(DTops[, c("chr","pos", newnames), with=F])
        }else{
          DT0 <- data.table(chr=chr, pos=pos)
          for(col in Vcols){
            DT0[DTops, (col):=get(col), on=c("chr","pos")]
          }
          return(DT0)
        }
    }
    
    if(closegds) gdsfmt::closefn.gds(gds)
    invisible(gc())
    data.table()
}

