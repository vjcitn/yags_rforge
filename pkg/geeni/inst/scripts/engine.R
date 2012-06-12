chunkable = function(data, n=NULL, chunksize=10) {
  # don't use a cursor but define an index:
  # stateless query resolution
  require(bit)
  if (is.null(n)) n= nrow(data)
  inds = chunk(from=1, to=n, by=chunksize, maxindex=n)
  datafun <- function(which) {
    if (which < 1 | which > length(inds)) return(NULL)
    if (which > length(inds)) stop("request for nonexistent chunk")
    tin = function(x) as.integer(as.hi(x))
    data[ tin(inds[[which]]), ]
    }
}

.dfbychunk = setRefClass("dfbychunk", fields = list(chunksize="integer", nchunk="integer",
    data="data.frame"),
  methods = list(show=function(...) cat(nrow(data), "x", ncol(data), "dfbychunk instance, chunksize", chunksize, "\n"),
  getchunk = function(which) chunkable(data, chunksize=chunksize)(which),
  initialize = function(..., chunksize, data) {
# do we really need this?
   if (missing(chunksize)) chunksize <<- 10L
   else chunksize <<- chunksize
   if (!missing(data)) data <<- data
   nchunk <<- as.integer(ceiling(nrow(data)/chunksize))
  }
))

combi = function (x, y)
   lapply(1:length(x), function(i) x[[i]] + y[[i]])

parols <- function (fmla, dfbyc, family, binit, maxit = 20, tol = 1e-06) 
{
    beta = binit
    del = Inf
    curit = 1
    nchunk = dfbyc$nchunk
    while (max(abs(del)) > tol) {
        delcomp = foreach(i = 1:nchunk, .combine = combi) %dopar% {
            fm = model.frame(fmla, dfbyc$getchunk(i))
            xi = model.matrix(fmla, fm)
            yi = model.response(fm)
	    resi = yi - xi%*%beta
            list(xpx = t(xi)%*%xi, xpr = t(xi)%*%resi, ssr=sum(resi^2), ni=length(yi))
            }
        xpxinv = solve(delcomp[[1]])
        del = xpxinv %*% delcomp[[2]]
        beta = beta + del
        curit = curit + 1
        if (curit > maxit) 
            stop(paste("maxit [", maxit, "] iterations exceeded"))
    }
    ssr = delcomp[[3]]
    N = delcomp[[4]]
    se = sqrt(diag((ssr/(N-length(del)))*xpxinv))
    list(stats=cbind(beta,se), ssr=ssr, N=N)
}
