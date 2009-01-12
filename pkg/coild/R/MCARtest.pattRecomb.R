
# 
# the idea is that the data are regimented enough
# to support resampling or permutation by cluster without
# 
# outcome: matrix nsubj x maxni
# 
# covariate: list of p matrices nsubj x maxni, one
#   per covariate
# 
# times: matrix nsubj x maxni of observation times
# 
# stratum: matrix nsubj x maxni of stratum indices
# 
# if we have a way of splitting and fitting by stratum,
# then permuting the rows of the stratum matrix
# gets what we want -- this would randomly reassign
# individuals to different strata, preserving number
# of individuals per stratum.
# 
# to proceed, convert from a standard thin gee dataset
# to the structure described here, where there is a 
# conforming vector of stratum indices
# 

demoSet <- function() {
y <- rnorm(60)
x1 <- runif(60)
x2 <- rnorm(60)
id <- sort(rep(1:15,4))
id[5] <- 1
id[9] <- 2
id[53:60] <- 59
sid <- split(id,id)
times <- unlist(lapply(sid,function(x)1:length(x)))
stratum <- sort(rep(1:3,20))
thin <- data.frame(y,x1,x2,id,times,stratum)
thin
}

thin2mats <- function(thin,xnames,yname="y",idname="id", timesname="times",
   stratumname="stratum") {
   pad <- function(x,n)
     {
     if (length(x)>n) stop("input longer than target length")
     if (length(x)==n) return(x)
     else { np <- n-length(x); return(c(x,rep(NA,np))) }
     }
   nx <- length(xnames)
   id <- thin[[idname]]
   lens <- unlist(lapply(split(id,id),length))
   nsub <- length(lens)
   maxni <- max(lens)
   ymat <- matrix(NA,nr=nsub,nc=maxni)
   tmat <- matrix(NA,nr=nsub,nc=maxni)
   idmat <- matrix(NA,nr=nsub,nc=maxni)
   svec <- rep(NA,nsub)
   y <- thin[[yname]]
   times <- thin[[timesname]]
   strat <- thin[[stratumname]]
   sy <- split(y,id)
   si <- split(id,id)
   st <- split(times,id)
   ss <- split(strat,id)
# process y, times, strat
   for (i in 1:nsub)
     {
     ymat[i,] <- pad(sy[[i]],maxni)
     tmat[i,] <- pad(st[[i]],maxni)
     idmat[i,] <- pad(si[[i]],maxni)
     svec[i] <- ss[[i]][1]
     }
# process x
   xin <- list()
   xout <- list()
   for (i in 1:nx)
     {
     tmp <- split( thin[[ xnames[i] ]], id )
     tmpmat <- matrix(NA, nr=nsub, nc=maxni)
     for (j in 1:nsub)
        tmpmat[j, ] <- pad( tmp[[j]], maxni)
     xout[[ xnames[i] ]] <- tmpmat
     }
   Key <- c(ymat=yname, tmat=timesname, idmat=idname, svec=stratumname)
   tmp <- list(ymat=ymat, xout=xout, tmat=tmat, idmat=idmat, stratvec=svec, Key=Key)
   names(tmp) <- c(yname,"xout",timesname,idname,stratumname,"Key")
   tmp
   }
 
assemByStrat <- function( ldslist )
{
 Z <- ldslist
 Key <- Z$Key
 ymat <- Z[[ Key["ymat"] ]]
 tmat <- Z[[ Key["tmat"] ]]
 idmat <- Z[[ Key["idmat"] ]]
 svec <- Z[[ Key["svec"] ]]
 nsubj <- nrow(ymat)
 nstrat <- length(us <- unique(svec))
 sustra <- sort(us)
 nobs <- sum(!is.na(ymat))
 ncthin <- 1+length(Z$xout)+1+1+1
 xnames <- names(Z$xout)
 nx <- length(xnames)
 znam <- names(Z)
 ynam <- znam[1]; tnam <- znam[3]; idnam <- znam[4]; snam <- znam[5]
 out <- matrix(NA,nr=nobs,nc=ncthin)
 dimnames(out) <- list(NULL,c(ynam,xnames,tnam,idnam,snam))
 curOutRow <- 1
 na.rm <- function(x) x[!is.na(x)]
 for (j in 1:nstrat)
  {
  rows2get <- (1:nsubj)[ svec==sustra[j] ]
  cury <- ymat[rows2get,,drop=FALSE]
  curid <- idmat[rows2get,,drop=FALSE]
  curt <- tmat[rows2get,,drop=FALSE]
  newy <- na.rm(t(cury))
  newid <- na.rm(t(curid))
  newt <- na.rm(t(curt))
  nnew <- length(newy)
  last <- curOutRow+nnew-1
  targs <- curOutRow:last
  out[targs,ynam] <- newy
  out[targs,idnam] <- newid
  out[targs,tnam] <- newt
  out[targs,snam] <- rep(sustra[j],nnew)
  for (k in 1:nx)
   {
   tmp <- Z$xout[[k]]
   curx <- tmp[rows2get,,drop=FALSE]
   newx <- na.rm(t(curx))
   out[targs,xnames[k]] <- newx
   }
  curOutRow <- last+1
  }
 data.frame(out)
}

splitGee <- function (formula, idname = "id", corstr, family, stratDf, stratName = "stratum")
{
    require(gee)
    if (!is.data.frame(stratDf))
        stop("need data frame in stratDf")
    ndat <- splitDf(stratDf, stratDf[[stratName]])
    nspl <- length(ndat)
    ans <- list()
    for (i in 1:nspl) {
        assign("tmpid", ndat[[i]][[idname]], pos = 1)
        assign("tmpdat", ndat[[i]], pos = 1)
        on.exit({
            rm(tmpid, pos = 1)
            rm(tmpdat, pos = 1)
        })
        ans[[i]] <- gee(formula, id = tmpid, corstr = corstr,
            family = family, data = tmpdat)
    }
    ans
}

 
# just a demonstration!

geeCoSize <- function (geefit,full) 
{
    B <- matrix(geefit$coef[5:6], nc = 1) - matrix(full$coef[5:6], nc=1)
    V <- geefit$robust.variance[5:6,5:6] + full$robust.variance[5:6,5:6]
    t(B) %*% ginv(V) %*% B
}

# more general implementation

coSizeMaker <- function (inds)
 {
  function (geefit, full) 
  {
     B <- matrix(geefit$coef[inds], nc = 1) - matrix(full$coef[inds], 
         nc = 1)
     V <- geefit$robust.variance[inds, inds] + full$robust.variance[inds, 
         inds]
     (t(B) %*% ginv(V) %*% B)
  }
 }


splitDf <- function( df, fac )
 {
 inds <- 1:nrow(df)
 sinds <- split(inds,fac)
 lens <- unlist(lapply(sinds,length))
 out <- list()
 for (i in 1:length(lens))
  out[[i]] <- df[sinds[[i]],]
 out
 }

#permByMissPatt <- function( df, misspatt )
#{
# Nrec <- nrow(df)
# if (length(misspat)!=Nrec) stop("misspatt and data frame nonconforming")
# nuniqmpat <- length(u <- sort(unique(misspatt)))
# if (any(u != 1:nuniqmpat)) stop("misspat values should be in 1:length(uniq(misspatt))
# per <- sample(1:Nrec,size=Nrec,replace=FALSE)
# df[per,]
#}
 
rperm <- function(x) sample(x,size=length(x),replace=FALSE)

MCARtest.pattRecomb <- function (formula, idname = "id", 
    family = gaussian, corstr = "independence",
    stratDf, yname = "y", xnames = c("x1", "x2"), timesname = "times",
    stratName = "stratum", sizeFun = geeCoSize, B = 20)
{
#
# this does not really implement the Chen Little method directly
# but follows it by measuring variation in stratum-specific
# estimators vs collapsed estimator
#
# the reference distribution is over permutations of stratum
# labels on clusters
#
    require(gee)
    fullid <<- stratDf[[idname]]
    on.exit(rm(fullid, pos = 1))
    full <- gee(formula = formula, id = fullid, family = family,
        corstr = corstr, data = stratDf)
    OBS <- splitGee(formula = formula, idname = idname, corstr = corstr,
        family = family, stratDf = stratDf, stratName = stratName)
    nstr <- length(OBS)
    ans.obs <- 0
    for (i in 1:nstr) ans.obs <- ans.obs + sizeFun(OBS[[i]],
        full)
    M <- thin2mats(stratDf, xnames = xnames, yname = yname, timesname = timesname,
        idname = idname, stratumname = stratName)
    ans.perm <- rep(0, B)
    for (i in 1:B) {
        Mtmp <- M
        Mtmp[[stratName]] <- rperm(Mtmp[[stratName]])
        newdf <- assemByStrat(Mtmp)
        pans <- try(splitGee(formula = formula, idname = idname,
            corstr = corstr, family = family, stratDf = newdf,
            stratName = stratName))
        if (!inherits(pans, "try-error")) {
            for (j in 1:nstr) ans.perm[i] <- ans.perm[i] + sizeFun(pans[[j]],
                full)
        }
    }
    list(ans.obs = ans.obs, ans.perm = ans.perm)
}


