\name{p.isRandSamp.MC}
\alias{p.isRandSamp.MC}
\title{
test that one vector is a random draw from another
}
\description{ obtain p-value of a montecarlo test for the hypothesis
that elements of sample x are a random draw from a population X}
\usage{
p.isRandSamp.MC(x, X, stat=mean, NSIM=19, inseed=123, 
  extreme=c("low", "high")[1])
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ vector representing the putative sample }
  \item{X}{ vector representing the putative population }
  \item{stat}{ statistic used to evaluate sampling hypothesis }
  \item{NSIM}{ number of iterations for montecarlo testing }
  \item{inseed}{ PRNG seed }
  \item{extreme}{ "low" or "high" }
}
\details{
}
\value{
%two-sided p-value
}
\references{ DLZ 11.2 }
\author{ VJ Carey }
\examples{
x <- rnorm(20)
X <- rnorm(200)
p.isRandSamp.MC( x, X )
}
\keyword{ sampling }
\keyword{ hypothesis testing }
\keyword{models}

