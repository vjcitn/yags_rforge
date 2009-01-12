\name{lagsAreCommon}
\alias{lagsAreCommon}
\title{ predicate on sampling times in a longitudinal dataset }
\description{
is true if intervisit gaps are the same for all subjects
in a longitudinal dataset
}
\usage{
lagsAreCommon(mat)
}
\arguments{
  \item{mat}{ an instance of class longimat }
}
\details{
}
\value{
boolean
}
\author{ VJ Carey }

\examples{
data(demod)
data(demod.times)
demoLD <- makeLongimat( mat=demod, times=demod.times )
lagsAreCommon(demoLD@times)
}
\keyword{ longitudinal data }
\keyword{ missing data }
\keyword{models}

