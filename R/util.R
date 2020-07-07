.. <- function (x, env=parent.frame()) 
{
  stopifnot(inherits(x, "character"))
  stopifnot(length(x) == 1)
  get(x, envir = parent.env(env))
}


set_colclass <- function(x, funcs){
  for(i in names(funcs)){
   if(i %in% colnames(x)){
     f <- funcs[[i]]
     if(!is.function(f) && is.character(f)){
       f <- get(f)
     }
     x[, (..("i")):=..("f")(get(..("i")))]
   } 
  }
}

submatrix <- function(x, rows=rownames(x), cols=colnames(x), missing=NA){
  stopifnot(!is.null(rownames(x)))
  stopifnot(!is.null(colnames(x)))
  
  res <- matrix(missing, nrow=length(rows), ncol=length(cols),
                dimnames = list(rows, cols))
  rowsinx <- intersect(rows, rownames(x))
  colsinx <- intersect(cols, colnames(x))
  
  if(length(rowsinx)>0 & length(colsinx) > 0){
    res[rowsinx, colsinx] <- as.matrix(x[rowsinx, colsinx])
  }
  
  res
}


strsplit2 <- function (str, split, fixed, perl, regex, coll, charclass, ...)
{
  type <- c(regex = !missing(regex), fixed = !missing(fixed), perl = !missing(perl),
            coll = !missing(coll), charclass = !missing(charclass))
  
  if(sum(type) == 0)
    type <- c(regex=TRUE)
  
  if (sum(type) != 1)
    stop("you have to specify either `regex`, `fixed`, `coll`, `perl`, or `charclass`")
  
  fun <- list(
    perl=stringi::stri_split_regex,
    regex=stringi::stri_split_regex,
    fixed=stringi::stri_split_fixed,
    coll=stringi::stri_split_coll,
    charclass=stringi::stri_split_charclass
  )[[names(type)[which(type)]]]
  
  fun(str, split, ...)
}

tstrsplit2 <- function (x, ..., fill = NA, type.convert = FALSE, keep, names = FALSE)
{
  ans = data.table::transpose(strsplit2(as.character(x), ...), fill = fill,
                              ignore.empty = FALSE)
  if (!missing(keep)) {
    keep = suppressWarnings(as.integer(keep))
    chk = min(keep) >= min(1L, length(ans)) & max(keep) <=
      length(ans)
    if (!isTRUE(chk))
      stop("'keep' should contain integer values between ",
           min(1L, length(ans)), " and ", length(ans), ".")
    ans = ans[keep]
  }
  if (type.convert)
    ans = lapply(ans, type.convert, as.is = TRUE)
  if (identical(names, FALSE))
    return(ans)
  else if (isTRUE(names))
    names = paste0("V", seq_along(ans))
  if (!is.character(names))
    stop("'names' must be TRUE/FALSE or a character vector.")
  if (length(names) != length(ans)) {
    str = if (missing(keep))
      "ans"
    else "keep"
    stop("length(names) (= ", length(names), ") is not equal to length(",
         str, ") (= ", length(ans), ").")
  }
  setattr(ans, "names", names)
  ans
}

#' paste cols of a data.table using key-value style
#' @param sep1 separate between key and value
#' @param sep2 separate between cols
#'
kvpaste <- function(df, sep1, sep2, cols){
  stopifnot(!missing(sep1))
  stopifnot(!missing(sep2))
  stopifnot(!missing(df))
  if (inherits(df, "matrix")) df <- as.data.table(df)
  #stopifnot(inherits(df, "data.frame"))
  if(missing(cols)) cols <- colnames(df)
  #df <- df[, lapply(.SD, as.character), .SDcols=cols]
  kvpaste_rcpp(df, sep1, sep2, nakey=FALSE)
}

#' split key-value strings into dataframe
#' @param x input strings
#' @param sep1 separate between key and value
#' @param sep2 separate between cols
#' @param na bool, FALSE for "", TRUE for NA_character_
kvsplit <- function(x, sep1, sep2, na=FALSE){
  out <- kvsplit_rcpp(x, sep1, sep2, na)
  setDT(out)
  out
}

split_long <- function(x, col, split, nm=col, keep.rowid=F, ...){
  stopifnot(is.data.table(x))
  stopifnot(all(col %in% colnames(x)))
  stopifnot(length(col) == length(nm))
  if(length(split)==1) split <- rep(split, length(col))
  stopifnot(length(col) == length(split))
  
  if(".rowid" %in% colnames(x)){
    stop("remove .rowid col first")
  }
  old.colorder <- colnames(x)
  
  x[[".rowid"]] <- 1L:nrow(x)
  setkey(x, .rowid)
  
  dt <- data.table()
  for(coli in col){
    s <- strsplit2(x[[coli]], split, ...)
    dti <- data.table(rep(1:length(s), sapply(s, length)), unlist(s, use.names=F))
    setnames(dti, c(".rowid", coli))
    if(nrow(dt)){
      dt <- unique(merge(dt, dti, all=T, allow.cartesian=T, by=".rowid"))
    }else{
      dt <- dti
    }
  }
  setnames(dt, c(".rowid", nm))
  if(setequal(c(".rowid", col), colnames(x))){
    x <- dt
  }else{
    x <- x[, setdiff(colnames(x), col), with = F]
    x <- x[dt, on=".rowid"]
  }
  
  if(!keep.rowid) {
    x <- x[, -".rowid", with=F]
    setcolorder(x, old.colorder)
  }else{
    setcolorder(x, c(".rowid", old.colorder))
  }
  x
}

na.replace <- function(x, y) {
  replace(x, is.na(x), y)
}

new.data.table <- function(fill, nrow, col.names) {
  dt <- data.table()
  if(length(col.names) > 0) {
    dt <- data.table(matrix(fill, nrow=nrow, ncol=length(col.names)))
    setnames(dt, col.names)
  }
  dt
}

