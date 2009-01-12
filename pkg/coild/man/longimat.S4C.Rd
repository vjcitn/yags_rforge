\name{longimat.S4C}
\alias{longimat.S4C}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ S4 class for longitudinal data, represented as an nxK matrix for
n subjects seen at K times }
\description{ S4 class for longitudinal data, represented as an nxK matrix for
n subjects seen at K times }
\usage{
#new( "longimat", ... )
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{id}{ numerical vector of subject identifiers }
  \item{labels}{ variable labels }
  \item{data}{ numerical response matrix }
  \item{times}{  matrix representing response times }
  \item{lagsAreCommon}{ logical indicating common intervisit
gaps for all subjects }
  \item{incomplete}{ logical indicating existence of some missing
responses }
  \item{monotone}{ logical indicating whether or not
missingness is monotone }
}
\details{
}
\value{
}
\references{ }
\author{ }
\note{ }

\seealso{  }

\examples{
data(demod)
data(demod.times)
demoLD <- makeLongimat( mat=demod, times=demod.times )
diggle.test(demoLD, function(x)x, stat=mean, NSIM=19, inseed=123)
}
\keyword{ longitudinal data }% at least one, from doc/KEYWORDS
\keyword{ missing data  }% __ONLY ONE__ keyword per line

