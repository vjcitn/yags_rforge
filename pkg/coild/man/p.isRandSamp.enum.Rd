\name{p.isRandSamp.enum}
\alias{p.isRandSamp.enum}
\title{ enumerative significance test of randomness of a sample }
\description{
exact significance level for Ho:x is a random sample from X
}
\usage{
p.isRandSamp.enum(x, X, stat=mean, extreme=c("low", "high")[1])
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ sample }
  \item{X}{ universe }
  \item{stat}{ statistic to be evaluated on x }
  \item{extreme}{ if "low", obtaining percentile of
observed stat in lower tail; otherwise 1- percentile in
upper tail}
}
\details{
}
\value{
}
\references{ Diggle 1989 Bcs }
\author{ VJ Carey }

\examples{
set.seed(1234)
X <- rnorm(25)
x <- sample(X,size=12)
# takes a long time!
# p.isRandSamp.enum(x,X)
# p.isRandSamp.enum(x,X,stat=median)
}
\keyword{ inference }
\keyword{models}

