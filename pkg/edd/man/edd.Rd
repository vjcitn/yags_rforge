\name{edd}
\alias{edd}
\title{ new expression density diagnostics interface }
\description{
this will replace edd.unsupervised; has more
sensible parameters
}
\usage{
edd(eset, distList=eddDistList, tx=c(sort,flatQQNormY)[[1]],
	refDist=c("multiSim", "theoretical")[1], 
	method=c("knn", "nnet", "test")[1], nRowPerCand=100, ...)
}
\arguments{
  \item{eset}{ eset -- instance of Biobase
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} class} 
  \item{distList}{ distList -- list comprised of eddDist objects}
  \item{tx}{ tx -- transformation of data and reference prior
to classification }
  \item{refDist}{ refDist -- type of reference distribution system to use}
  \item{method}{ method -- type of classifier to use.  knn is k-nearest
neighbors, nnet is neural net, test is max p-value from ks.test}
  \item{nRowPerCand}{ nRowPerCand -- number of realizations for
a multiSim reference system}
  \item{\dots}{ \dots -- parameters to classifiers}
}
\details{
Classifies genes according to distributional shape, by
comparing observed expression distributions to
a collection of references, which may be simulated
or evaluated theoretically.

The distList argument is important.  It enumerates the
catalog of distributions for classification of gene
expression vectors by distributional shape.  See the
HOWTO-edd vignette for information on how this list
is constructed and how it can be extended.

The tx argument specifies how the data are processed
for comparison to the reference catalog.  This is
a function on a vector returning a vector, but the
input and the output need not have the same length.
The default value of tx is sort, which entails that
the order statistics are treated as multivariate data
for classification.  

The refDist argument selects the type of reference
catalog.  Options are 'multiSim', for which the
reference consists of nRowPerCand realizations of
each catalog entry, and 'theoretical', for which
the reference consists of one vector of quantiles
for each catalog entry.

The method argument selects the type of classifier.
It would be desirable to allow this to be a function,
but there is insufficient structure on classifier
argument and return value structure to permit this
at present; see the e1071 package for some work
on handling various classifiers programmatically
(e.g., \code{tune}).

}
\value{
a character vector or factor depending on the classifier
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


\seealso{\code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}  }

\examples{
require(Biobase)
data(sample.ExpressionSet)
# should filter to genes with reasonable variation
table( edd(sample.ExpressionSet, meth="nnet", size=10, decay=.2) )
library(golubEsets)
data(Golub_Merge)
madvec <- apply(exprs(Golub_Merge),1,mad)
minvec <- apply(exprs(Golub_Merge),1,min)
keep <- (madvec > median(madvec)) & (minvec > 300)
gmfilt <- Golub_Merge[keep==TRUE,]
ALL <- gmfilt$ALL.AML=="ALL"
gall <- gmfilt[,ALL==TRUE]
gaml <- gmfilt[,ALL==FALSE]
alldists <- edd(gall, meth="nnet", size=10, decay=.2)
amldists <- edd(gaml, meth="nnet", size=10, decay=.2)
table(alldists,amldists)
amldists2 <- edd(gaml, meth="nnet", refDist="theoretical", size=10, decay=.2)
table(amldists,amldists2)
}
\keyword{ models }
