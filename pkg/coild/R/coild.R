
.initClasses <- function(env) {

setClass("dataset", representation( id="numeric", labels="character", data="matrix"),
 where=env)

setClass("longimat", contains="dataset", representation(times="matrix", 
			lagsAreCommon="logical", incomplete="logical",
    			monotone="logical"),where=env)

setClass("ildConf", representation(n="numeric",
   K="numeric", subjectModelName="character", subjectModelParms="list",
   subjectNonrespFun="function", monotone="logical", 
   subjectRealizer="function"), where=env)

setClass("ildRealization", contains="ildConf", representation=
    representation(mat="longimat", seed="numeric"), where=env)

#setClass("annotateL", representation( name="character",
#         fname="character", key="numeric"), where=env)


#assign("LlU95", new("annotateL", name="LlU95", fname="u95LL", key=2),
#                 env=env)

setGeneric("is.monomiss", function(object) standardGeneric("is.monomiss"),
   where=env)

setMethod("is.monomiss", "longimat", function(object) { 
 x <- object@data
 is.monomiss.mat(x)
 },where=env)


setGeneric("getni",  function(object) standardGeneric("getni"), where=env)

# following needs more to prevent getting a dataset
# with missing baseline (e.g., min(j)>1)
#setMethod("[", c("longimat","numeric"), function(x, i ,
#j ,..., drop=TRUE)
# {
# if (missing(j)&missing(i)) return(x)
# if (missing(j)) 
#       return(makeLongimat(x@data[i,,drop=drop], x@times[i,,drop=drop]))
# if (missing(i)) 
#       return(makeLongimat(x@data[,j,drop=drop], x@times[,j,drop=drop]))
# return(makeLongimat(x@data[i,j,drop=drop], x@times[i,j,drop=drop]))
# })

setMethod("getni", "longimat",
 function(object) {
  x <- object@data
  apply(x,1,function(x)sum(!is.na(x)))
 }, where=env)

#setGeneric("dim", where=env)

setMethod("dim", "longimat", function(x) dim(x@data), where=env)

setGeneric("toGEEformat", function(object) standardGeneric("toGEEformat"), where=env)

setMethod("toGEEformat", "longimat", function(object) {
 x <- object
 nsub <- nrow(x@data)
 if (length(x@id) != nsub) stop("id and data elements nonconforming")
 ntim <- ncol(x@data)
 id <- rep(x@id,rep(ntim,nsub))
 y <- as.numeric(t(x@data))
 tim <- as.numeric(t(x@times))
 data.frame(y=y,tim=tim,id=id)
 }, where=env)

setGeneric("sampleFrom", function(object)standardGeneric("sampleFrom"), where=env)

setMethod("sampleFrom", "ildConf", function(object) {
 saveSeed <- .Random.seed
 mat <- matrix(NA, nr=object@n, nc=object@K)
 times <- mat
 for (i in 1:object@n) 
   {
   mat[i,] <- evalSlot(object, "subjectRealizer")
   times[i,] <- seq(1,object@K)
   }
 new("ildRealization", mat=makeLongimat(mat, times),
    K=object@K, n=object@n,
    subjectModelName=object@subjectModelName, subjectModelParms=
    object@subjectModelParms,
   subjectNonrespFun=object@subjectNonrespFun, monotone=object@monotone,
    seed=saveSeed, subjectRealizer=object@subjectRealizer)
 }, where=env)

setGeneric("makeMonotone",  function(object) standardGeneric("makeMonotone"), where=env)

setMethod("makeMonotone", "longimat", function(object) {
 nr <- nrow(object@data)
 z <- object@data
 for (i in 1:nr)
  z[i,] <- make.monotone.vec(object@data[i,])
 return(makeLongimat(z,object@times))
 }, where=env)

} # end of .initClasses

make.monotone.vec <- function(x) { 
if (!any(is.na(x))) return(x)
n <- length(x); 
first.NA <- min((1:n)[is.na(x)]);
   return (c(x[1:(first.NA-1)],rep(NA,n-first.NA+1)))
}


runif(1)

lmisspatt <- function (x) 
{
    logi <- !is.na(x)
    char <- logi
    char <- matrix(as.character(char), nc = ncol(logi))
    char[char == "TRUE"] <- "p"
    char[char == "FALSE"] <- " "
    tmp <- apply(char, 1, function(x) paste(x, sep = "", collapse = ""))
    t(t(table(tmp)))
}

make.monotone.vec <- function(x) { 
if (!any(is.na(x))) return(x)
n <- length(x); 
first.NA <- min((1:n)[is.na(x)]);
   return (c(x[1:(first.NA-1)],rep(NA,n-first.NA+1)))
}

is.monomiss.vec <- function(vec) {
 if (all(!is.na(vec))) return(TRUE)
 if (all(is.na(vec))) return(TRUE)
 n <- length(vec)
 first.miss <- min((1:n)[is.na(vec)])
 last.nonmiss <- max((1:n)[!is.na(vec)])
 if (last.nonmiss > first.miss) return(FALSE)
 return(TRUE)
 }

is.monomiss.mat <- function(mat) {
 x <- mat
 all(apply(x,1,is.monomiss.vec))
 }

lagsAreCommon <- function(mat) 
 {
 nr <- nrow(mat)
 del1 <- diff(mat[1,])
 eqd <- function(x,y)all(diff(x)==y)
 all(apply(mat,1,eqd,del1))
 }

makeLongimat <- function(mat,times,labels=paste("V",1:ncol(mat),sep=""),id=1:nrow(mat)) {
 if (missing(times)) stop("times must be a matrix of observation times\nconforming to mat")
 if (!is.matrix(mat)) stop("mat must be a matrix")
 if (!all(dim(mat) == dim(times))) stop("mat and times do not conform")
 incomplete <- any(is.na(mat))
 monotone <- is.monomiss.mat(mat)
 lagsAreCommon <- lagsAreCommon(times)
 new("longimat", data=mat, labels=labels, times=times, id=id,
   incomplete=incomplete, monotone=monotone, lagsAreCommon=lagsAreCommon)
 }

ind.Rr <- function(obj, rel=">=") {
 rel <- rel
 ni <- getni(obj)
 function(k) {
  x <- obj@data
  (1:nrow(x))[do.call(rel,list(ni, k))]
 }
}

# IND.R <- ind.Rr(demo, ">=")
# IND.r <- ind.Rr(demo, "==")

p.isRandSamp.enum <- function (x, X, stat = mean, extreme = c("low", "high")[1]) 
{
    m <- stat(x)
    N <- length(X)
    n <- length(x)
    nch <- nCm(length(X), length(x))
    tmp <- rep(NA, nch)
    allsel <- combn(X, length(x))
    for (i in 1:nch) {
        tmp[i] <- stat(allsel[, i])
    }
    if (extreme == "low") 
        n.past <- sum(m < tmp)
    else n.past <- sum(m >= tmp)
    n.past/nch
}

p.isRandSamp.MC <- function(x, X, stat=mean, NSIM=19, 
    inseed=123, extreme=c("low", "high")[1]) {
# MC test for hypothesis that elements of x are a random
# sample from the population X
 set.seed(inseed)
 m <- stat(x)
 N <- length(X)
 n <- length(x)
 tmp <- rep(NA,NSIM)
 for (i in 1:NSIM)
  {
  tmp[i] <- stat(sample(X,size=n))
  }
 if (extreme == "low") 
     n.past <- sum(m < tmp)
 else
     n.past <- sum(m >= tmp)
 prop <- n.past/(NSIM+1)
# if (prop > .5) return(2*(1-prop))
# else return(prop)
 prop
}

last.only <- function(x) x[length(x)]
 
diggle.test <- function(obj, scorefn, stat=mean, NSIM=19, inseed=123,
  method=c("MC","enum")[1], extreme="low", require.monotone=TRUE) {
 if (require.monotone & !obj@monotone) stop("must have monotone missingness")
 if (!obj@lagsAreCommon) stop("must have lagsAreCommon visit times")
 if (!any(is.na(obj@data))) stop("no missing data")
 ind.R <- ind.Rr(obj)
 ind.r <- ind.Rr(obj,"==")
 Y <- obj@data
 K <- ncol(Y)
 tmp <- rep(NA,K-1)
 for (i in 1:(K-1))
  {
  Rind <- ind.R(i)
  rind <- ind.r(i)
  if (length(rind)==length(Rind)) next
  if (length(rind)==0) next # no missing at that visit
  H <- apply(Y[Rind,1:i,drop=FALSE], 1, scorefn)
  h <- apply(Y[rind,1:i,drop=FALSE], 1, scorefn)
  if (method == "enum")
     tmp[i] <- p.isRandSamp.enum( h, H, stat=stat, extreme=extreme)
  else
     tmp[i] <- p.isRandSamp.MC( h, H, stat=stat, NSIM=NSIM, inseed=inseed, extreme=extreme )
  }
 tmp
}

evalSlot <- function(x, slot, ...)
 slot(x, slot)(x, ...)

demoSubjectRealizer <- function(this, ...)
 {
 tmp <- do.call( this@subjectModelName, this@subjectModelParms )
 this@subjectNonrespFun( tmp )
 }

demoNRfun <- function(x,ub=1.5)
 {
 n <- length(x)
 xx <- (x - median(x))/mad(x)
 if (!any(z <- xx[2:n]>ub)) return(x)
 drop <- min((1:n)[c(FALSE,z)])
 x[(drop):n] <- NA
 x
 }

makeC1Reg <- function(a,b,k,sig)
 {
 x <- 1:k
 a+b*x+rnorm(k,0,sig)
 }

#C1parms <- list(a=0,b=1,k=5,sig=.5)

#C1 <- new("ildConf", n=100, K=5, subjectModelName="makeC1Reg",
# subjectModelParms = C1parms, subjectNonrespFun = demoNRfun,
# subjectRealizer=demoSubjectRealizer)
#SC1 <- sampleFrom(C1)

WuFollman4 <- function (x) 
{
# Biometrics March 1999 p77
    require(nlme)
    df <- toGEEformat(x)
    m1 <- lme(fixed = y ~ tim, random = ~tim | id, data = df, na.action=na.omit)
    ree <- m1$coef$random$id
    X <- 1*is.na(df$y)
    bi1 <- as.numeric(ree[,1])
    bi2 <- as.numeric(ree[,2])
    TT <- df$tim  # should have no missingness
    newdf <- data.frame(x=X, tim=TT, bi1=bi1, bi1t=TT*bi1, bi2=bi2, bi2t=TT*bi2,id=x@id)
    gee(x~TT+bi1+bi1t+bi2+bi2t,fam=binomial, data=newdf, id=id)
}

