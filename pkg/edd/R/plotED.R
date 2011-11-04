plotED <- function (x, y, data = NULL, is.centered = FALSE, 
    ...) 
{
    dname <- paste("d", stub(x), sep = "")
    argl <- list()
    ps <- parms(x)
    nps <- names(ps)
    for (i in 1:length(ps <- parms(x))) {
        argl[[nps[i]]] <- ps[i]
    }
    dfun <- function(x) {
        argl[["x"]] <- x
        do.call(dname, argl)
    }
    curve(dfun, from = plotlim(x)[1], to = plotlim(x)[2], ylab = "density", 
     xlab = ifelse(is.null(data),"x","centered/scaled data relocated to nominal support"), 
	main=tag(x), ...)
    if (!is.null(data)) {
        if (is.centered) 
            tdata <- Mad(x) * data + med(x)
        else {
            tdata <- Mad(x) * centerScale(data) + med(x)
        }
        lines(density(tdata), lty = 2)
    }
}
