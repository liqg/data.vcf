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
     x[, (get("i")):=..("f")(get(..("i")))]
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

