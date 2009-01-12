\name{makeLongimat}
\alias{makeLongimat}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
generator function for instances of class \code{longimat}
}
\usage{
makeLongimat(mat, times, labels=paste("V", 1:ncol(mat), sep = ""), id)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{mat}{ matrix of responses with one row per subject and
one column per repeated observation }
  \item{times}{ a conforming matrix of observation times }
  \item{labels}{ column labels }
  \item{id}{ subject labels }
}
\details{
generates elements from a longitudinal data representation particularly
convenient for assessment of the associated missing data process
}
\value{
an object of S4 class \code{longimat}
}
\author{ VJ Carey }

\examples{
data(demod)
data(demod.times)
demoLD <- makeLongimat( mat=demod, times=demod.times )
}
\keyword{ longitudinal data }
\keyword{ missing data }
\keyword{models}

