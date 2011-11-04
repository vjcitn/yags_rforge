latEDtable <- function(x,baselist=eddDistList,reorder=NULL) {
 ltags <- unlist(lapply(baselist,latexTag))
 ctags <- unlist(lapply(baselist,tag))
 names(ltags) <- ctags
 names(ctags) <- ltags
 if(length(dd <- dim(x))==2) {
   xn <- dimnames(x)
   names(xn) <- NULL
   mx <- matrix(x,nr=dd[1],nc=dd[2])
   dimnames(mx) <- lapply(xn,function(x)ltags[x])
 }
 else
  {
   xn <- names(x)
   mx <- matrix(x,nr=1)
   dimnames(mx) <- list(" ",ltags[xn])
  }
 if (is.null(reorder)) return(mx)
 else if (length(dd)==2) return(mx[reorder,reorder])
 else return(mx[,reorder,drop=FALSE])
}
  
