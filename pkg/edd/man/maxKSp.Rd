\name{maxKSp}
\alias{maxKSp}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ classify densities by maximum KS p-value }
\description{
compares a sample to reference distributions, identifying
the closest reference by maximum KS p
}
\usage{
maxKSp( x, is.centered=TRUE, dists=eddDistList, thresh=.1 )
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ matrix of samples, distributions to be classified by row }
  \item{is.centered}{ have data been brought to median 0, mad 1 }
  \item{dists}{ list of instances of class eddDist }
  \item{thresh}{ p-value above which some test
must lie, or else 'outlier' is declared }
}
\details{
}
\value{
a vector of classifications
}

\examples{
X <- rbind( rnorm(30), runif(30), rchisq(30,1))
tX <- t(apply(X,1,centerScale))
apply(tX,1,maxKSp)
}
\keyword{ models }% at least one, from doc/KEYWORDS
