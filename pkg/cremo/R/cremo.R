#
# source code for cremo package
#

#
# S4 phase: define methods that extract model information
# from R fitted model objects
# the ultimate purpose: to have a set of methods that
# obtain coefficient values, variable names, s.e.s, p-values,
# and extra goodness of fit information.  all these functions
# return named numeric vectors, with the exception of getVarn
#
setGeneric("getCoef", function(x) standardGeneric("getCoef"))
setGeneric("getVarn", function(x) standardGeneric("getVarn"))
setGeneric("getSE", function(x) standardGeneric("getSE"))
setGeneric("getP", function(x) standardGeneric("getP"))
setGeneric("getExtra", function(x) standardGeneric("getExtra"))

setMethod("getCoef", "lm", function(x) coef(x))
setMethod("getCoef", "geese", function(x) x$beta)
setMethod("getCoef", "gee", function(x) x$coefficient)
setMethod("getCoef", "alr", function(x) x$coefficient)
setMethod("getCoef", "yagsResult", function(x) 
          yags.glmReport(x)[,1])

setMethod("getSE", "lm", function(x) summary(x)$coef[,2] )
setMethod("getSE", "geese", function(x) summary(x)$mean[,2] )
setMethod("getSE", "gee", function(x) summary(x)$coef[,4] )
setMethod("getSE", "alr", function(x) summary(x)$coef[,2] )
setMethod("getSE", "yagsResult", function(x) 
          yags.glmReport(x)[,4])

setMethod("getP", "lm", function(x) summary(x)$coef[,4] )
setMethod("getP", "geese", function(x) summary(x)$mean[,4] )
setMethod("getP", "gee", function(x) 2*pnorm(-abs(summary(x)$coef[,5] )))
setMethod("getP", "alr", function(x) 2*pnorm(-abs(summary(x)$coef[,3] )))
setMethod("getP", "yagsResult", function(x) 2*pnorm(-abs(yags.glmReport(x)[,5] )))

setMethod("getVarn", "lm", function(x) names(coef(x)))
setMethod("getVarn", "geese", function(x) names(x$beta))
setMethod("getVarn", "gee", function(x) names(coef(x)))
setMethod("getVarn", "alr", function(x) names(coef(x)))
setMethod("getVarn", "yagsResult", function(x) row.names(yags.glmReport(x)))

setMethod("getCoef", "clogit", function(x) exp(coef(x)))
setMethod("getVarn", "clogit", function(x) names(coef(x)))
setMethod("getSE", "clogit", function(x) rep(NA, length(coef(x))))
setMethod("getP", "clogit", function(x) summary(x)$coef[,5])
setMethod("getExtra", "clogit", function(x) summary(x)$sctest )

setMethod("getExtra", "lm", function(x) {
  sx <- summary(x)
  c(r.sq=sx$r.sq, se.resid=sx$sigma) })
setMethod("getExtra", "geese", function(x) {
  sx <- summary(x)
  sxc <- as.numeric(unlist(sx$corr))
  c(cor=sxc[1], se.cor=sxc[2]) })
setMethod("getExtra", "gee", function(x) {
  sx <- summary(x)
  c(cor=sx$working.corr[1,2], se.cor=NA) })
setMethod("getExtra", "alr", function(x) {
  sx <- summary(x)
#
# here we need a regimented naming convention
# and to interleave estimates and s.es
  avec <- sx$alpha[,1]
  sevec <- sx$alpha[,2]
  nparms <- length(avec)
  ans <- as.numeric(rbind(avec,sevec))
  names(ans) <- paste(c("alp", "sealp"), rep(1:nparms,each=2),
    sep=".")
  ans
})
setMethod("getExtra", "yagsResult", function(x) {
  tmp <- yags.adeqReport(x)
  ans <- c(substr(x@corstruct.tag,1,4), round(c(x@alpha[1], rj1=tmp@rj1, rj2=tmp@rj2),3) )
  names(ans) <- c("mod", "alp[1]", "RJ1", "RJ2")
  ans
  })
  


concatMods <- function(x, mstub="mod",type=c("se","p")[1],dig=3,clean=TRUE) {
#
# the workhorse of cremo -- takes a list of model objects and
# produces a character matrix with columns corresponding to
# models and rows corresponding to parameters or extra material
#
  nx = names(x)
  padsubvec <- function (subv, full, pad="-") 
  {
#
# pad a (names) vector subv that only partly matches full with
# value of pad in each non-matched element ... probably
# could be done better 
#
    minds <- match(subv, full)
    ok <- (1:length(full)) %in% minds
    tmp <- full
    tmp[ok == FALSE] <- pad
    tmp
  }
  selWpad <- function(x,y,pad="-")
  {
#
# suppose names vector y is padded by padsubvec; the elements
# in x that correspond to non-padded names are kept, and the
# remainder are set to value of pad
#
    tx <- c(x,padval=NA)
    names(tx)[length(tx)] <- pad 
    tx[y]
  }
  interleave <- function(x,y) {
#interleave rows of character matrices
    if(!(nrow(x)==nrow(y) & ncol(x)==ncol(y)))stop("x y must have same dims")
    if (!is.character(x[1,1])) stop("x must be character")
    if (!is.character(y[1,1])) stop("y must be character")
    nr <- 2*nrow(x)
    out <- matrix(" ", nr=nr, nc=ncol(x))
    j <- 1
    for (i in seq(1,nr,2))
     {
     out[i,] <- x[j,]
     out[i+1,] <- y[j,]
     j <- j+1
     }
    rn <- rownames(x)
    rn <- rbind(rn," ")
    dimnames(out) <- list(as.character(rn),dimnames(x)[[2]])
    out
   }
  clnup <- function (cm) 
   {
#
# get rid of unnecessary NA printings
#
       tm <- t(matrix(as.character(t(cm)), nr = ncol(cm)))
       tm[is.na(tm)] <- " "
	 gg = grep("NA", tm)
        if (length(gg)>0) tm[gg] <- " "
       tm[tm == "(p=NA)"] <- " "
       tm[tm == "(SE=NA)"] <- " "
       dimnames(tm) <- dimnames(cm)
       tm
   }
  pwrap <- function(x) paste("(p=",x,")",sep="")
  sewrap <- function(x) paste("(SE=",x,")",sep="")
#
# now the processing of the model objects begins
#
  coList <- lapply(x, getCoef)
  seList <- lapply(x, getSE)
  PList <- lapply(x, getP)
  eList <- lapply(x, getExtra)
  varnList <- lapply(x, getVarn)
  nm <- length(coList)
  for (i in 1:nm)  # number of models
   {
   names(coList[[i]]) <- varnList[[i]]
   names(seList[[i]]) <- varnList[[i]]
   names(PList[[i]]) <- varnList[[i]]
   names(varnList[[i]]) <- varnList[[i]]
   }
allvar <- unique(unlist(varnList))
allEx <- unique(unlist(en <- lapply(eList,names)))
#
# out will hold the names of coefficients properly padded
# to agree with indices in the vector of the union 
# of all variable names (allvar), which is then used
# to populate the {se,p,co}mat
#
out <- matrix(0, nr=length(allvar), nc=nm)
semat <- pmat <- comat <- matrix(NA, nr=length(allvar), nc=nm)
if (is.null(nx)) modn <- paste(mstub,1:nm,sep="")
else modn = nx
eout <- emat <- matrix(NA, nr=length(allEx), nc=nm)
#
# eout will hold names of all the extra material, and
# emat will hold the selected results
#
for (i in 1:nm)
  {
  out[,i] <- padsubvec( varnList[[i]], allvar )
  comat[,i] <- selWpad(coList[[i]],out[,i])
  pmat[,i] <- selWpad(PList[[i]],out[,i])
  semat[,i] <- selWpad(seList[[i]],out[,i])
  eout[,i] <- padsubvec(en[[i]],allEx)
  emat[,i] <- selWpad(eList[[i]],eout[,i])
  }
dimnames(comat) <- dimnames(pmat) <- list(allvar, modn)
cpmat <- pmat
csemat <- semat
#cpmat[] <- pwrap(as.character(round(pmat,dig)))
cpmat[] <- pwrap(formatC(pmat,dig=dig))
csemat[] <- sewrap(as.character(round(semat,dig)))
comat[] <- as.character(round(comat,dig))
if (is.numeric(emat)) emat[] <- as.character(round(emat,dig))
else emat[] <- as.character(emat)
rownames(emat) <- allEx #names(eList[[1]])
if (type == "p") ans <- interleave(comat,cpmat)
else if (type == "se") ans <- interleave(comat,csemat)
ans <- rbind(ans, "--", emat)
dans <- data.frame(ans)
#
# these are 'unique' unobtrusive row.names for blank rows
#
bl <- c("", " ", "  ", "   ", "    ","     ", "      ", "        ",
"         ",".", ". "," .", ".  ")
bl <- c(bl, paste("-", bl, sep=""))
rn <- rownames(ans)
bln <- rn %in% c("", " ")
rn[bln] <- bl[1:sum(bln)]
rownames(dans) <- rn
if (!clean) return(dans)
else return(clnup(dans))
}

