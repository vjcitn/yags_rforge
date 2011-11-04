\name{plotED}
\alias{plotED}
\title{visualize an eddDist object  }
\description{
plots an eddDist model; can superimpose density sketch
from transformed data
}
\usage{
plotED(x, y, data=NULL, is.centered=FALSE, ...)
}
\arguments{
  \item{x}{ eddDist object}
  \item{y}{ not used}
  \item{data}{ optional vector of data to be
superimposed in the form of a kernel density estimate}
  \item{is.centered}{ is.centered: logical indicating
that the data vector has zero median and unit mad}
  \item{\dots}{ \dots: options passed to curve}
}
\author{Vince Carey <stvjc@channing.harvard.edu>}



\examples{
#
# show the first 8 supplied reference dists
if (interactive()){
omf <- par()$mfrow
oas <- par()$ask
on.exit(par(mfrow=omf,ask=oas))
par(mfrow=c(4,2))
par(ask=TRUE)
}
set.seed(1234)
for (i in 1:8) plotED(eddDistList[[i]])
# illustrate the superposition
if (interactive()) par(mfrow=c(1,1))
x <- rnorm(30,3,4)
plotED(N01,data=x) # relocates/scales x
y <- 12*rbeta(30,2,8)+4
plotED(B28,data=y)
}
\keyword{ models }
