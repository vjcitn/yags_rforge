\name{ind.Rr}
\alias{ind.Rr}
\title{ closure generator, supporting discrimination of
imminent dropouts and persistent subjects }
\description{
The isolation of subsets of persistent subjects and imminent
dropouts at each visit time is a bit cumbersome.  Let obj be a
longimat isntance.  Indices
of subjects who persist to visit k+1 <= K are in the
vector ind.Rr(obj, rel=">=")(k), and indices of those
who drop out after visit k but before k+1 are in
the vector ind.Rr(obj, rel="==")(k).
}
\usage{
ind.Rr(obj, rel=">=")
}
\arguments{
  \item{obj}{ an instance of the longimat class }
  \item{rel}{ either ">=" or "==" }
}
\details{
The aim is to have a function with domain k = 1,...,K-1
and range $R_k$ (indices of subjects present beyond visit k)
or range $r_k$ (indices of subjects dropping out just
after visit k).  ind.Rr is a closure generator yielding
functions, e.g., ind.R(k) or ind.r(k) to be used during
iteration over k.
}
\value{
returns functions to be used during iteration over the argument;
see the example.
}
\author{ VJ Carey }
\examples{
data(demod)
data(demod.times)
demoLD <- makeLongimat( mat=demod, times=demod.times )
ind.R <- ind.Rr(demoLD)
ind.r <- ind.Rr(demoLD, rel="==")

ind.R(2)
ind.r(2) # should be a subset of ind.R(2)
}

\keyword{ longitudinal data }
\keyword{ incomplete data }
\keyword{models}

