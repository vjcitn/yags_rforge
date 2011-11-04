
pmixnorm <- function(x,p1,m1,s1,m2,s2)
   p1 * pnorm(x,m1,s1) + (1-p1) * pnorm(x, m2, s2)
qmixnorm <- function(pct,p1,m1,s1,m2,s2)
   sapply(pct, function(z) optimize(function(x)
     (pmixnorm(x,p1,m1,s1,m2,s2)-z)^2,c(-20,20))$minimum)
rmixnorm <- function(n,p1,m1,s1,m2,s2)
     qmixnorm(runif(n),p1,m1,s1,m2,s2)
dmixnorm <- function (x, p1, m1, s1, m2, s2)
{
    if (p1 >= 1)
        stop("mixture fraction must be less than 1")
    p1 * dnorm(x, m1, s1) + (1 - p1) * dnorm(x, m2, s2)
}


makeCandmat.theor <- function (nPerRow, distList, pctiles=
  ppoints(nPerRow))
{
#
# each row is the sorted rescaledQuantiles of the
# 9 reference distributions N01, chisq1 t3 ln01, u01, mix1, mix2, beta2,8, beta8,2
#
    rescaledQuantiles <- function(eddist, pctiles) {
        (qfun(eddist)(pctiles) - med(eddist))/Mad(eddist)
    }
    t(sapply(distList, rescaledQuantiles, pctiles ))
}

eddObsolete <- function(eset, 
    ref=c("multiCand","uniCand","test","nnet")[1], k=10, l=6,
    nnsize=6, nniter=200)
 {
#
# oddities of this function
#  1) it uses the flatQQnorm transformation for data when
#     it should suffice to just use the edf
#  2) it has odd defaults and some hardcoded options inside
#  3) the logic is strange
#
 require(class)
 require(nnet)
 nsamp <- ncol(exprs(eset))
 if (ref=="multiCand" || ref=="nnet") 
      {
      refmat <- makeCandmat.raw(nPerRow=nsamp,nRowPerCand=100)
      refmat <- t(apply(refmat,1,centerScale))
      fq.ref <- fq.matrows(refmat)
      if (ref == "nnet")
        {
        Y <- row.names(fq.ref)
        row.names(fq.ref) <- NULL ##so we don't get a warning next
        dfcx <- data.frame(Y=Y,fq.ref)
        refnet <- nnet(Y~.,data=dfcx,size=nnsize,iter=nniter, decay=0.01)
        }
      }
 else if (ref=="uniCand")
      fq.ref <- makeCandmat.theor(nPerRow=nsamp)
 else if (ref=="test")
      return( t(apply(exprs(eset), 1, function(x)maxKSp(centerScale(x)))) )
 fq.test <- fq.matrows(t(apply(exprs(eset),1,centerScale)))
 if (ref == "nnet")
   return(predict(refnet,newdata=data.frame(fq.test),type="class"))
 else knn( train=fq.ref, cl=row.names(fq.ref), test=fq.test, k=k, l=l )
 }

edd <- function( eset, distList=eddDistList, tx=c(sort, flatQQNormY)[[1]],
  refDist=c("multiSim", "theoretical")[1],
  method=c("knn", "nnet", "test")[1], nRowPerCand=100, ...) {
    nsamp <- ncol(exprs(eset))
    testSet <- data.frame(t(apply(exprs(eset),1,function(x)tx(centerScale(x)))))
    if (method=="test") return(t(apply(testSet,1,maxKSp, dists=distList, ...)))
    if (refDist == "multiSim")
      ref <- makeCandmat.raw(nPerRow=nsamp, nRowPerCand=nRowPerCand,
                             dists=distList)
    else
      ref <- makeCandmat.theor(nPerRow=nsamp, distList=distList)
    ref <- t(apply(ref,1,tx))
    if (method == "knn")
      ans <- knn( train=ref, cl=row.names(ref), 
		test=testSet, ...)
    else
        {
        Y <- row.names(ref)
	dimnames(ref) <- list(NULL,NULL)
        dfcx <- data.frame(Y=Y,ref)
        names(dfcx) <- c("Y",names(testSet))
        refnet <- nnet(Y~.,data=dfcx,...)
	ans <- predict(refnet, newdata=testSet, type="class")
        }
    ans
}
	
			  

mkt <- function() {
# a test matrix with samples from some of the reference dists
test <- matrix(NA,nr=140,nc=50)
test[1:20,] <- rnorm(1000)
test[21:40,] <- rchisq(1000,1)
test[41:60,] <- rt(1000,3)
test[61:80,] <- rlnorm(1000)
test[81:100,] <- runif(1000)
test[101:120,] <- rmixnorm(1000,.75,0,1,4,1)
test[121:140,] <- rbeta(1000,2,8)
test
}

# labels for the test matrix entries
testcl <- c(rep("n01",20), rep("csq1",20), rep("t3",20), 
      rep("ln",20), rep("u",20), rep("mix1",20),rep("b28",20))

"fq.matrows" <-
function(mat)
 t(apply(mat,1,function(x,...){tmp <- flatQQNorm(x); tmp$y[order(tmp$x)]}))

#"ctr" <-
#function(x) (x-median(x))/mad(x)

"flatQQNorm" <-
function (y)
{
    n <- length(y)
    x <- qnorm(ppoints(n))[order(order(y))]
    y <- y - x
    list(x = x, y = y)
}
"flatQQNormY" <- function(x) flatQQNorm(x)$y

#"matbox.strat" <-
#function (x, strat, ...) 
#{
## split a matrix by strat, then plot boxplots
## for each strat (separate boxplots within panel for each stratum)
#    spar <- par()
#    p <- ncol(x)
#    m <- list()
#    cl <- strat
#    cltags <- unique(cl)
#    ncl <- length(cltags)
#    for (i in 1:ncl) {
#        m[[i]] <- list()
#        for (j in 2:p) m[[i]][[j - 1]] <- x[cl == cltags[i], 
#            j]
#        n <- length(m[[i]][[1]])
#        boxplot(m[[i]], main = paste("class", cltags[i], "size=", 
#            n), ...)
#    }
#    par(spar)
#    invisible(0)
#}

centerScale <- function(x) (x-median(x))/mad(x)

rmixnorm.alt <- function(n,p1,m1,s1,m2,s2) 
  c( rnorm(floor(p1*n),m1,s1), rnorm(ceiling((1-p1)*n),m2,s2) )

pmixnorm <- function(x,p1,m1,s1,m2,s2)  {
  if (p1 >= 1) stop("mixture fraction must be less than 1")
  p1*pnorm(x,m1,s1)+(1-p1)*pnorm(x,m2,s2)
  }

#
# maxKSp is a good illustration of what should be  done
# generically to compare a vector to a collection of
# distributional models
#
# in this setting we have 3 possible comparisons
# a and b treat edf(x) as an n-dimensional multivariate datum
# and use multivariate classification methods like knn or nnet
#  a) compare edf(x) to theoretical edf 
#  b) compare edf(x) to edfs of simulated 
#            realizations from theoretical dists
#  c) compare edf(x) to theoretical edf using KS test
#
# problem: to specify conveniently and flexibly the
# various theoretical edfs of interest so that the user
# can add new ones or omit some in the default collection
#
# for maxKSp we need to be able to do do.call("ks.test", list
# (x, pfun, parms) ) and we may have to rescale x to match
# the median and mad for pfun
#
# for makeCandmat.raw we need to do do.call(gen,parmlist)
#
# for makeCandmat.theor we need to evaluate quantiles at
# selected probabilities
#
# thus an eddDist object would have a character stub
# for the name (e.g., "norm" or "t") from which the
# generator, quantile function or cdf can be derived
# by pasting, a numeric parameter vector assumed to
# be conformable for the distribution in hand, numerical
# precomputed values of median and mad
#
setClass("eddDist", representation(stub="character",
		parms="numeric", median="numeric", mad="numeric",
		tag="character", plotlim="numeric", latexTag="character"))
setGeneric("stub", function(x)standardGeneric("stub"))
setMethod("stub", "eddDist", function(x)x@stub)
setGeneric("plotlim", function(x)standardGeneric("plotlim"))
setMethod("plotlim", "eddDist", function(x)x@plotlim)
setGeneric("tag", function(x)standardGeneric("tag"))
setMethod("tag", "eddDist", function(x)x@tag)
setGeneric("latexTag", function(x)standardGeneric("latexTag"))
setMethod("latexTag", "eddDist", function(x)x@latexTag)
setGeneric("parms", function(x)standardGeneric("parms"))
setMethod("parms", "eddDist", function(x)x@parms)
setGeneric("med", function(x)standardGeneric("med"))
setMethod("med", "eddDist", function(x)x@median)
setGeneric("Mad", function(x)standardGeneric("Mad"))
setMethod("Mad", "eddDist", function(x)x@mad)
setGeneric("CDFname", function(x)standardGeneric("CDFname"))
setMethod("CDFname", "eddDist", function(x) paste("p",stub(x),sep=""))
setGeneric("qfName", function(x)standardGeneric("qfName"))
setMethod("qfName", "eddDist", function(x) paste("q",stub(x),sep=""))
setGeneric("qfun", function(e)standardGeneric("qfun"))
setMethod("qfun", "eddDist", function(e)  {
  argl <- list()
  pn <- names(parms(e))
  for (i in 1:length(parms(e)))
     argl[[pn[i]]] <- parms(e)[i]
# following assumes arg to quantile func is always named p
  function(x) { argl[["p"]] <- x; do.call(qfName(e),argl) }
})
setGeneric("genName", function(x)standardGeneric("genName"))
setMethod("genName", "eddDist", function(x) paste("r",stub(x),sep=""))
setGeneric("testVec", function(x,eddd,is.centered)standardGeneric("testVec"))
setMethod("testVec", c("numeric","eddDist","logical"), 
  function(x,eddd,is.centered) {
    if (is.centered) x <- Mad(eddd)*x + med(eddd)
    argl <- list(x, CDFname(eddd))
    if (length(Parms <- parms(eddd))>0)
       for (j in 1:length(Parms)) argl[[2+j]] <- Parms[j]
 oldopt = options()
 options(warn=-1)
 on.exit(options(oldopt))
    do.call(ks.test, argl)
})

N01 <- new("eddDist", stub="norm", parms=c(mean=0,sd=1), 
  median=0, mad=1, tag="N(0,1)", plotlim=c(-3,3),
  latexTag="$\\Phi$")
CS1 <- new("eddDist", stub="chisq", parms=c(df=1), median=.46, mad=.6,
 tag="X^2(1)", plotlim=c(0.001,6), latexTag="$\\chi^2_1$")
T3 <- new("eddDist", stub="t", parms=c(df=3), median=.0, mad=1.15,
 tag="t(3)", plotlim=c(-4,4), latexTag="$t_3$")
LN01 <- new("eddDist", stub="lnorm", parms=c(meanlog=0,sdlog=1), 
   median=1, mad=.88, tag="logN(0,1)", plotlim=c(0,5),latexTag="$LN_{0,1}$")
U01 <- new("eddDist", stub="unif", parms=c(min=0,max=1), median=.5, mad=.37,
 tag="U(0,1)", plotlim=c(-.1,1.1), latexTag="$U_{0,1}$")
MIXN1 <- new("eddDist", stub="mixnorm", parms=c(p1=.75,m1=0,s1=1,m2=4,s2=1),
     median=.42, mad=1.57, tag=".75N(0,1)+.25N(4,1)", plotlim=c(-3,7),
	latexTag="$\\frac{3}{4}\\Phi+\\frac{1}{4}\\Phi_{4,1}$")
MIXN2 <- new("eddDist", stub="mixnorm", parms=c(p1=.25,m1=0,s1=1,
     m2=4,s2=1),
     median=3.58, mad=1.57, tag=".25N(0,1)+.75N(4,1)", plotlim=c(-3,7),
	latexTag="$\\frac{1}{4}\\Phi+\\frac{3}{4}\\Phi_{4,1}$")
B28 <- new("eddDist", stub="beta", parms=c(shape1=2,shape2=8), 
	median=.18, mad=.12, tag="B(2,8)", plotlim=c(0,1),
	latexTag="$\\beta_{2,8}$")
B82 <- new("eddDist", stub="beta", parms=c(shape1=8,shape2=2), 
	median=.82, mad=.12, tag="B(8,2)", plotlim=c(0,1),
	latexTag="$\\beta_{8,2}$")

# put these in order of similarity
eddDistList <- list(N01=N01, T3=T3, LN01=LN01, CS1=CS1, 
   B82=B82, U01=U01, B28=B28, MIXN1=MIXN1, MIXN2=MIXN2)

maxKSp <- function( x, is.centered=TRUE, dists=eddDistList, thresh=.1 ) {
 ps <- sapply(dists, function(z) testVec(x,z,is.centered)$p.value)
 if (all(ps<thresh)) return("outlier") # need to handle doubt as well
 ind <- (1:length(ps))[ps>=max(ps)][1]
 return(tag(dists[[ind]]))
}


makeCandmat.raw <- function(nPerRow=20, nRowPerCand=20, dists=
  eddDistList, centerScale=TRUE)
{
 ncand <- length(dists)
 out <- matrix(NA, nr=nRowPerCand*ncand, nc=nPerRow)
 gnames <- rep(NA,ncand)
 for (i in 1:ncand)
  {
  gnames[i] <- tag(dists[[i]])
  inds <- ((i-1)*nRowPerCand+1):(i*nRowPerCand)
  n2gen <- nRowPerCand*nPerRow
  argl <- list(n2gen)
  for (j in 1:length(pl <- parms(dists[[i]])))
      argl[[1+j]] <- pl[j]
  out[inds,] <- do.call(genName(dists[[i]]), argl)
  if (centerScale) {
    me <- med(dists[[i]])
    ma <- Mad(dists[[i]])
    out[inds,] <- (out[inds,]-me)/ma
    }
  }
 row.names(out) <- rep(gnames,each=nRowPerCand)
 out
}


