".." <- function (x) 
{
  stopifnot(inherits(x, "character"))
  stopifnot(length(x) == 1)
  get(x, parent.frame(4))
}

set_cols <- function(x, funcs){
  for(i in names(funcs)){
   if(i %in% colnames(x)){
     f <- funcs[[i]]
     if(!is.function(f)){
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