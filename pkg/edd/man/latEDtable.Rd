\name{latEDtable}
\alias{latEDtable}
\title{ use latex tags for dimnames of table created from edd classification }
\description{use latex tags for dimnames of table created from edd classification
}
\usage{
latEDtable(x, baselist=eddDistList, reorder=NULL)
}
\arguments{
  \item{x}{ x -- a table (possibly two dimensional) of results
of applying edd}
  \item{baselist}{ baselist -- a list of eddDist objects}
  \item{reorder}{ reorder -- a numeric vector describing
how to re order the table rows/columns}
}
\details{
for use with xtable rendering.  table will give
results with margin names in lexicographic order;
reorder can be used to get a different order.
}
\value{
a matrix with dimnames computed from the latexTag slots
of the corresponding eddDist objects
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
require(Biobase)
data(sample.ExpressionSet)
# should filter to genes with reasonable variation
rawTab <- table( edd(sample.ExpressionSet, meth="nnet", size=10, decay=.2) )
latEDtable(rawTab)
library(xtable)
xtable(latEDtable(rawTab))
#
realTags <- sapply(eddDistList,tag)
reo <- match(realTags,names(rawTab))
xtable(latEDtable(rawTab,reorder=reo))
}
\keyword{ models }
